#!/usr/bin/env python

"""
    epitopepredict analysis methods for workflows
    Created September 2013
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, shutil, string, types
import csv, glob, pickle, itertools
import math
import re
import time, random
from collections import OrderedDict
from operator import itemgetter
import numpy as np
import pandas as pd
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from . import base, sequtils, tepitope, utilities, peptutils

home = os.path.expanduser("~")
#fix paths!
genomespath = os.path.join(home, 'epitopedata')
datadir = os.path.join(home, 'testpredictions')

def get_AAcontent(df, colname, amino_acids=None):
    """Amino acid composition for dataframe with sequences"""
    return df.apply( lambda r: peptutils.get_AAfraction(str(r[colname]), amino_acids), 1)

def net_charge(df, colname):
    """Net peptide charge for dataframe with sequences"""
    return df.apply( lambda r: peptutils.net_charge(r[colname]),1)

def isoelectric_point(df):
    def getpi(seq):
        X = ProteinAnalysis(seq)
        return X.isoelectric_point()
    return df.apply( lambda r: getpi(r.peptide),1)

def peptide_properties(df, colname='peptide'):
    """Find hydrophobicity and net charge for peptides"""

    df['hydro'] = get_AAcontent(df, colname)
    df['net_charge'] = net_charge(df, colname)
    return df

def _center_nmer(x, n):
    """Get n-mer sequence for a peptide centered in the middle.
    This should be applied to a dataframe per row.
    Returns: a single sequence centred on the peptide"""

    seq = x['translation']
    size = x.end-x.start
    l = int((size-n)/2.0)
    if size>n:
        if size%2 == 1: l1 = l+1
        else: l1=l
        start = x.start+l1
        end = x.end-l
    elif size<=n:
        if size%2 == 1: l1 = l-1
        else: l1=l
        start = x.start+l1
        end = x.end-l
    if start<=0:
        d=1-start
        start = start+d
        end = end+d
    seq = seq[start:end]
    #print(size, x.peptide, x.start, x.end, l, l1, start, end, seq, len(seq))
    return seq

def _split_nmer(x, n, key, margin=3):
    """Row based method to split a peptide in to multiple n-mers
        if it's too large. Returns a dataframe of 3 cols so should be
        applied using iterrows and then use concat"""

    size = x.end-x.start
    m = margin

    if size <= n+m:
        seq = _center_nmer(x, n)
        return pd.DataFrame({'peptide': seq,
                             'start':x.start,'end':x.end},index=[0])
    else:
        seq = x[key]
        o=size%n
        #print (x.peptide, size, o)
        if o<=margin:
            size=size-o
            seq = _center_nmer(x, size)
            print (size)
        seqs=[]
        S=[];E=[]
        if x.start==0: s=1
        else: s=0
        for i in range(s, size, n):
            if i+n>size:
                seqs.append(seq[o:o+n])
                S.append(x.start+o)
                E.append(x.start+o+n)
            else:
                seqs.append(seq[i:i+n])
                S.append(x.start+i)
                E.append(x.start+i+n)
        seqs = pd.Series(seqs)
        d = pd.DataFrame({'peptide':seqs,'start':S,'end':E})
        return d

def get_nmer(df, genome, length=20, seqkey='peptide', how='center', margin=3):
    """
    Get n-mer peptide surrounding a set of sequences using the host
    protein sequence.
    Args:
        df: input dataframe with sequences
        genome: genome dataframe with host sequences
        length: length of nmer to return
        seqkey: column name of sequence to be processed
        how: method to create the n-mer, 'split' will try to split up
            the sequence into overlapping n-mes of length is larger than size
        margin: allow
    Returns:
        pandas Series with nmer values
    """

    cols = ['locus_tag','gene','translation','length']
    cols = list(set(cols) & set(genome.columns))
    #merge with genome dataframe but must keep index for re-merging!
    temp = df.merge(genome[cols],left_on='name',right_on='locus_tag',
                    how='left').set_index(df.index)
    if not 'end' in list(temp.columns):
        temp = base.get_coords(temp)

    temp = base.get_coords(temp)
    if how == 'center':
        res = temp.apply( lambda r: _center_nmer(r, length), 1)
    elif how == 'split':
        res=[]
        for n,r in temp.iterrows():
            d = _split_nmer(r, length, seqkey, margin)
            res.append(d)
            d['index']=n
            d.set_index('index',inplace=True)
        res = pd.concat(res)
    return res

def create_nmers(df, genome, key='nmer', length=20, margin=1):
    """
    Add n-mers to a dataframe of sequences by splitting them up
    and updating the start/end coords of each row in the dataframe
    """

    x = get_nmer(df, genome, how='split', length=length, margin=margin)
    x = x.rename(columns={'peptide':key})
    df = df.drop(['start','end'],1)
    x = df.merge(x,left_index=True,right_index=True).reset_index(drop=True)
    return x

def get_overlaps(df1, df2, label='overlap', how='inside'):
    """
    Overlaps for 2 sets of sequences where the positions in host sequence are stored
    in each dataframe as 'start' and 'end' columns
    Args:
        df1 : first set of sequences, a pandas dataframe with columns called
                start/end or pos
        df2: second set of sequences
        label: label for overlaps column
        how: may be 'any' or 'inside'
    Returns:
        First DataFrame with no. of overlaps stored in a new column
    """

    new=[]
    a = base.get_coords(df1)
    b = base.get_coords(df2)

    def overlap(x,y):
        f=0
        #print x['name'],x.peptide
        #print (x.start,x.end)
        for i,r in y.iterrows():
            if how == 'inside':
                if ((x.start<=r.start) & (x.end>=r.end)):
                    f+=1
            elif how == 'any':
                if ((x.start<r.start) & (x.end>r.start)) or \
                   ((x.start>r.start) & (x.start<r.end)):
                    #t = abs(r.start-x.start)
                    #print (a, b)
                    f+=1
            #print (r.start,r.end, f)
        return f

    for n,df in a.groupby('name'):
        found = b[b.name==n]
        df[label] = df.apply(lambda r: overlap(r,found),axis=1)
        new.append(df)
    result = pd.concat(new)
    print ('%s with overlapping sequences' %len(result[result[label]>0]))
    return result

def get_orthologs(seq, db=None, expect=1, hitlist_size=400, equery=None,
                  email=''):
    """
    Fetch orthologous sequences using remote or local blast and return the records
    as a dataframe.
    Args:
        seq: sequence to blast
        db: the name of a local blast db
        expect: expect value
        equery: Entrez Gene Advanced Search options,
                (see http://www.ncbi.nlm.nih.gov/books/NBK3837/)
    Returns:
        blast results in a pandas dataframe
    """

    from Bio.Blast import NCBIXML,NCBIWWW
    from Bio import Entrez, SeqIO
    Entrez.email = email

    print ('running blast..')
    if db != None:
        #local blast
        SeqIO.write(SeqRecord(Seq(seq)), 'tempseq.faa', "fasta")
        sequtils.local_blast(db, 'tempseq.faa', output='my_blast.xml', maxseqs=100)
        result_handle = open("my_blast.xml")
        df = sequtils.get_blast_results(result_handle)
    else:
        try:
            result_handle = NCBIWWW.qblast("blastp", "nr", seq, expect=expect,
                                  hitlist_size=500,entrez_query=equery)
            time.sleep(2)
            savefile = open("my_blast.xml", "w")
            savefile.write(result_handle.read())
            savefile.close()
            result_handle = open("my_blast.xml")
            df = sequtils.get_blast_results(result_handle, local=False)
        except Exception as e:
            print ('blast timeout')
            return

    df = df.drop(['subj','positive','query_length','score'],1)
    df.drop_duplicates(subset=['definition','perc_ident'], inplace=True)
    df = df[df['perc_ident']!=100]
    return df

def alignment_to_dataframe(aln):
    alnrows = [[str(a.id),str(a.seq)] for a in aln]
    df = pd.DataFrame(alnrows,columns=['accession','seq'])
    return df

def align_blast_results(df, aln=None, idkey='accession', productkey='definition'):
    """
    Get gapped alignment from blast results using muscle aligner.
    """

    sequtils.dataframe_to_fasta(df, idkey=idkey, seqkey='sequence',
                        descrkey=productkey, outfile='blast_found.faa')
    aln = sequtils.muscle_alignment("blast_found.faa")
    alnrows = [[a.id,str(a.seq)] for a in aln]
    alndf = pd.DataFrame(alnrows,columns=['accession','seq'])
    #res = df.merge(alndf, left_index=True, right_index=True)
    res = df.merge(alndf, on=['accession'])
    res = res.drop('sequence',1)
    #get rid of duplicate hits
    #res.drop_duplicates(subset=['definition','seq'], inplace=True)
    res = res.sort_values(by='identity',ascending=False)
    print ('%s hits, %s filtered' %(len(df), len(res)))
    return res, aln

def get_species_name(s):
    """Find [species name] in blast result definition"""

    m = re.search(r"[^[]*\[([^]]*)\]", s)
    if m == None:
        return s
    return m.groups()[0]

def find_conserved_sequences(seqs, alnrows):
    """
    Find if sub-sequences are conserved in given set of aligned sequences
    Args:
        seqs: a list of sequences to find
        alnrows: a dataframe of aligned protein sequences
    Returns:
        a pandas DataFrame of 1 or 0 values for each protein/search sequence
    """

    f=[]
    for i,a in alnrows.iterrows():
        sequence = a.seq
        found = [sequence.find(j) for j in seqs]
        f.append(found)
    for n in ['species','accession','name']:
        if n in alnrows.columns:
            ind = alnrows[n]
            break
    s = pd.DataFrame(f,columns=seqs,index=ind)
    s = s.replace(-1,np.nan)
    s[s>0] = 1
    res = s.count()
    return s

def epitope_conservation(peptides, alnrows=None, proteinseq=None, blastresult=None,
                         blastdb=None, perc_ident=50, equery='srcdb_refseq[Properties]'):
    """
    Find and visualise conserved peptides in a set of aligned sequences.
    Args:
        peptides: a list of peptides/epitopes
        alnrows: a dataframe of previously aligned sequences e.g. custom strains
        proteinseq: a sequence to blast and get an alignment for
        blastresult: a file of saved blast results in plain csv format
        equery: blast query string
    Returns:
        Matrix of 0 or 1 for conservation for each epitope/protein variant
    """

    import seaborn as sns
    sns.set_context("notebook", font_scale=1.4)

    if alnrows is None:
        if proteinseq == None:
            print ('protein sequence to blast or alignment required')
            return
        if blastresult == None or not os.path.exists(blastresult):
            blr = get_orthologs(proteinseq, equery=equery, blastdb=blastdb)
            if blr is None:
                return
            #if filename == None: filename = 'blast_%s.csv' %label
            blr.to_csv(blastresult)
        else:
            blr = pd.read_csv(blastresult, index_col=0)
        #blr = blr[blr.perc_ident>=perc_ident]
        alnrows, aln = align_blast_results(blr)
        #print (sequtils.formatAlignment(aln))

    if 'perc_ident' in alnrows.columns:
        alnrows = alnrows[alnrows.perc_ident>=perc_ident]
    if 'definition' in alnrows.columns:
        alnrows['species'] = alnrows.definition.apply(get_species_name)
    c = find_conserved_sequences(peptides, alnrows).T

    c = c.dropna(how='all')
    c = c.reindex_axis(c.sum(1).sort_values().index)
    if len(c) == 0:
        print ('no conserved epitopes in any sequence')
        return
    return c

def _region_query(P, eps, D):
	neighbour_pts = []
	for point in D:
		if abs(P - point)<eps:
			neighbour_pts.append(point)
	return neighbour_pts

def _expand_cluster(P, neighbour_pts, C, c_n, eps, min_pts, D, visited):

    flatten = lambda l: [i for sublist in l for i in sublist]
    C[c_n].append(P)
    for point in neighbour_pts:
        if point not in visited:
            visited.append(point)
            neighbour_pts_2 = _region_query(point, eps, D)
            if len(neighbour_pts_2) >= min_pts:
                neighbour_pts += neighbour_pts_2
        #print (point,C)
        if point not in flatten(C):
            C[c_n].append(point)

def _dbscan(D, eps=5, minsize=2):
    """
    1D intervals using dbscan. Density-Based Spatial clustering.
    Finds core samples of high density and expands clusters from them.
    """
    from numpy.random import rand
    noise = []
    visited = []
    C = []
    c_n = -1
    for point in D:
        visited.append(point)
        neighbour_pts = _region_query(point, eps, D)
        if len(neighbour_pts) < minsize:
            noise.append(point)
        else:
            C.append([])
            c_n+=1
            _expand_cluster(point, neighbour_pts, C, c_n,eps, minsize, D, visited)

    C = [i for i in C if len(i)>=minsize]
    #for cl in C:
    #    print (cl)
    return C

def dbscan(B=None, x=None, dist=7, minsize=4):
    """Use dbscan algorithm to cluster binder positions"""

    if B is not None:
        if len(B)==0:
            return
        x = sorted(B.pos.astype('int'))
    clusts = _dbscan(x, dist, minsize)

    '''from sklearn.cluster import DBSCAN
    X = np.array(list(zip(x,np.zeros(len(x)))), dtype=np.int)
    db = DBSCAN(eps=dist, min_samples=minsize)
    db.fit(X)
    labels = db.labels_
    n_clusters_ = len(set(labels))
    clusts=[]
    for k in range(n_clusters_):
        my_members = labels == k
        #print "cluster {0}: {1}".format(k, X[my_members, 0])
        if len(X[my_members, 0])>0:
            clusts.append(list(X[my_members, 0]))'''

    #print (clusts)
    return clusts

def find_clusters(binders, dist=None, min_binders=2, min_size=12, max_size=50,
                 genome=None, colname='peptide'):
    """
    Get clusters of binders for a set of binders.
    Args:
        binders: dataframe of binders
        dist: distance over which to apply clustering
        min_binders : minimum binders to be considered a cluster
        min_size: smallest cluster length to return
        max_size: largest cluster length to return
        colname: name for cluster sequence column
    Returns:
        a pandas Series with the new n-mers (may be longer than the initial dataframe
        if splitting)
    """

    C=[]
    grps = list(binders.groupby('name'))
    length = binders.head(1).peptide.str.len().max()
    #print (length)
    if dist == None:
        dist = length+1
        #print ('using dist for clusters: %s' %dist)
    for n,b in grps:
        if len(b)==0: continue
        clusts = dbscan(b,dist=dist,minsize=min_binders)
        if len(clusts) == 0:
            continue
        for c in clusts:
            gaps = [c[i]-c[i-1] for i in range(1,len(c))]
            C.append([n,min(c),max(c)+length,len(c)])

    if len(C)==0:
        print ('no clusters')
        return pd.DataFrame()
    x = pd.DataFrame(C,columns=['name','start','end','binders'])
    x['length'] = (x.end-x.start)
    x = x[x['length']>=min_size]
    x = x[x['length']<=max_size]

    #if genome data available merge to get peptide seq
    '''cols = ['locus_tag','translation']
    if 'gene' in genome.columns:
        cols.append('gene')
    if genome is not None:
        x = x.merge(genome[cols],
                    left_on='name',right_on='locus_tag')
        x[colname] = x.apply(lambda r: r.translation[r.start:r.end], 1)
        x = x.drop(['locus_tag','translation'],1)
        x = x.drop_duplicates(colname)
    x = x.sort_values(by=['binders'],ascending=False)
    x = x.reset_index(drop=True)'''
    print ('%s clusters found in %s proteins' %(len(x),len(x.groupby('name'))))
    print
    return x

def randomized_lists(df, n=94, seed=8, filename='peptide_lists'):
    """
    Return a randomized lists of sequences from a dataframe. Used for
    providing peptide lists for assaying etc.
    """

    #cols for excel file
    lcols=['label','peptide']
    #nb always use same seed
    np.random.seed(seed=seed)
    plist = df.reset_index(drop=True)
    plist = plist.reindex(np.random.permutation(plist.index))
    plist['label'] = plist['name']+'_'+plist.index.astype(str)
    plist.to_csv('%s.csv' %filename)

    writer = pd.ExcelWriter('%s.xls' %filename)
    #chunks = np.array_split(plist.index, groups)
    i=1
    chunks = plist.groupby(np.arange(len(plist)) // n)
    for g,c in chunks:
        c[lcols].to_excel(writer,'list'+str(i))
        i+=1
    #also save by method for easy reference
    #for i,g in plist.groupby('method'):
    #    g.sort_values(by='name')[lcols].to_excel(writer,'method'+str(i))
    plist.to_excel(writer,'original data', float_format='%.2f')
    writer.save()

'''def genome_analysis(datadir,label,gname,method):
    """this method should be made independent of web app paths etc"""

    path = os.path.join(datadir, '%s/%s/%s' %(label,gname,method))
    #path='test'
    gfile = os.path.join(genomespath,'%s.gb' %gname)
    g = sequtils.genbank2Dataframe(gfile, cds=True)
    b = getAllBinders(path, method=method, n=5)
    P = base.getPredictor(method)
    res = b.groupby('name').agg({P.scorekey:[np.mean,np.size,np.max]}).sort()
    res.columns = res.columns.get_level_values(1)
    res = res.merge(g[['locus_tag','length','gene','product','order']],
                            left_index=True,right_on='locus_tag')
    res['perc'] = res['size']/res.length*100
    res = res.sort('perc',ascending=False)

    top = b.groupby('peptide').agg({P.scorekey:np.mean,'allele':np.max,
                    'name': lambda x: x}).reset_index()
    top = top.sort(P.scorekey,ascending=P.rankascending)
    cl = findClusters(b, method, dist=9, minsize=3)
    if cl is not None:
        gc = cl.groupby('name').agg({'density':np.max})
        res = res.merge(gc,left_on='locus_tag',right_index=True)
    #print res[:10]
    return res'''

def tmhmm(fastafile=None, infile=None):
    """
    Get TMhmm predictions
    Args:
        fastafile: fasta input file to run
        infile: text file with tmhmm prediction output
    """

    if infile==None:
        tempfile = 'tmhmm_temp.txt'
        cmd = 'tmhmm %s > %s' %(fastafile,tempfile)
        infile = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    tmpred = pd.read_csv(infile, delim_whitespace=True, comment='#',
                      names=['locus_tag','v','status','start','end'])
    tmpred = tmpred.dropna()
    print ('tmhmm predictions for %s proteins' %len(tmpred.groupby('locus_tag')))
    lengths=[]
    for i,row in tmpred.iterrows():
        if row.status == 'TMhelix':
            lengths.append(row.end-row.start)
    #print np.mean(lengths), np.std(lengths)
    return tmpred

def signalP(infile=None,genome=None):
    """Get signal peptide predictions"""

    if genome != None:
        seqfile = Genome.genbank2Fasta(genome)
        tempfile = 'signalp_temp.txt'
        cmd = 'signalp -t gram+ -f short %s > %s' %(seqfile,tempfile)
        infile = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    sp = pd.read_csv(infile,delim_whitespace=True,comment='#',skiprows=2,
                      names=['locus_tag','Cmax','cpos','Ymax','ypos','Smax',
                            'spos','Smean','D','SP','Dmaxcut','net'])
    #print sp[sp.SP=='Y']
    return sp

def get_seqdepot(seq):
    """Fetch seqdepot annotation for sequence"""

    from epitopepredict import seqdepot
    reload(seqdepot)
    sd = seqdepot.new()
    aseqid = sd.aseqIdFromSequence(seq)
    try:
        result = sd.findOne(aseqid)
    except Exception as e:
        print (e)
        result=None
    return result

def prediction_coverage(expdata, binders, key='sequence', perc=50, verbose=False):
    """
    Determine hit rate of predictions in experimental data
    by finding how many top peptides are needed to cover % positives
    Args:
        expdata: dataframe of experimental data with peptide sequence and name column
        binders: dataframe of ranked binders created from predictor
        key: column name in expdata for sequence
    Returns:
        fraction of predicted binders required to find perc total response
    """

    def getcoverage(data, peptides, key):
        #get coverage for single sequence
        target = math.ceil(len(data)*perc/100.0)
        if verbose == True:
            print (len(data), target)
        #print data[key]
        #print peptides[peptides.isin(data[key])]
        found=[]
        count=0
        for p in peptides:
            for i,r in data.iterrows():
                #print p, r[key]
                if r[key] in found:
                    continue
                if r[key].find(p)!=-1 or p.find(r[key])!=-1:
                    found.append(r[key])
                    if verbose == True:
                        print (count, p, r[key])
                    continue
            count+=1
            if len(found) >= target:
                if verbose == True:
                    print (count, target)
                    print ('--------------')
                return count
        if verbose == True:
            print ('not all sequences found', count, target)
        return count

    total = 0
    for name, data in expdata.groupby('name'):
        peptides = binders[binders.name==name].peptide
        if len(peptides) == 0:
            continue
        if verbose == True: print (name)
        #print binders[binders.name==name][:5]
        c = getcoverage(data, peptides, key)
        total += c

    #print (total, total/float(len(binders))*100)
    return round(total/float(len(binders))*100,2)

def test_features():
    """test feature handling"""

    fname = os.path.join(datadir,'MTB-H37Rv.gb')
    df = sequtils.genbank2Dataframe(fname, cds=True)
    df = df.set_index('locus_tag')
    keys = df.index
    name='Rv0011c'
    row = df.ix[name]
    seq = row.translation
    prod = row['product']
    rec = SeqRecord(Seq(seq),id=name,description=prod)
    fastafmt = rec.format("fasta")
    print (fastafmt)
    print (row.to_dict())
    ind = keys.get_loc(name)
    previous = keys[ind-1]
    if ind<len(keys)-1:
        next = keys[ind+1]
    else:
        next=None
    return

def testrun(gname):

    method = 'tepitope'#'iedbmhc1'#'netmhciipan'
    path='test'
    gfile = os.path.join(genomespath,'%s.gb' %gname)
    df = sequtils.genbank2Dataframe(gfile, cds=True)
    #names = list(df.locus_tag[:1])
    names=['VP24']
    alleles1 = ["HLA-A*02:02", "HLA-A*11:01", "HLA-A*32:07", "HLA-B*15:17", "HLA-B*51:01",
              "HLA-C*04:01", "HLA-E*01:03"]
    alleles2 = ["HLA-DRB1*0101", "HLA-DRB1*0305", "HLA-DRB1*0812", "HLA-DRB1*1196", "HLA-DRB1*1346",
            "HLA-DRB1*1455", "HLA-DRB1*1457", "HLA-DRB1*1612", "HLA-DRB4*0107", "HLA-DRB5*0203"]
    P = base.getPredictor(method)
    P.iedbmethod='IEDB_recommended' #'netmhcpan'
    P.predictProteins(df,length=11,alleles=alleles2,names=names,
                        save=True, path=path)
    f = os.path.join('test', names[0]+'.mpk')
    df = pd.read_msgpack(f)
    P.data=df
    #b = P.get_binders(data=df)
    #print b[:20]
    base.getScoreDistributions(method, path)
    return

def test_conservation(label,gname):
    """Conservation analysis"""

    tag='VP24'
    pd.set_option('max_colwidth', 800)
    gfile = os.path.join(genomespath,'%s.gb' %gname)
    g = sequtils.genbank2Dataframe(gfile, cds=True)
    res = g[g['locus_tag']==tag]
    seq = res.translation.head(1).squeeze()
    print (seq)
    #alnrows = getOrthologs(seq)
    #alnrows.to_csv('blast_%s.csv' %tag)
    alnrows = pd.read_csv('blast_%s.csv' %tag,index_col=0)
    alnrows.drop_duplicates(subset=['accession'], inplace=True)
    alnrows = alnrows[alnrows['perc_ident']>=60]
    seqs=[SeqRecord(Seq(a.sequence),a.accession) for i,a in alnrows.iterrows()]
    print (seqs[:2])
    sequtils.distanceTree(seqs=seqs)#,ref=seqs[0])
    #sequtils.ETETree(seqs, ref, metric)
    #df = sequtils.getFastaProteins("blast_found.faa",idindex=3)
    '''method='tepitope'
    P = base.getPredictor(method)
    P.predictSequences(df,seqkey='sequence')
    b = P.get_binders()'''
    return

def find_conserved_peptide(peptide, recs):
    """Find sequences where a peptide is conserved"""

    f=[]
    for i,a in recs.iterrows():
        seq = a.sequence.replace('-','')
        found = seq.find(peptide)
        f.append(found)
    s = pd.DataFrame(f,columns=['found'],index=recs.accession)
    s = s.replace(-1,np.nan)
    #print s
    res = s.count()
    return s

def test():
    gname = 'ebolavirus'
    label = 'test'

    testrun(gname)
    #testBcell(gname)
    #testgenomeanalysis(label,gname,method)
    #testconservation(label,gname)
    #testFeatures()
    return

if __name__ == '__main__':
    pd.set_option('display.width', 600)
    test()
