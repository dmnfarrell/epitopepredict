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
#import matplotlib
#matplotlib.use('agg')
import pylab as plt
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

def getAAContent(df, amino_acids=None):
    """Amino acid composition for dataframe with sequences"""
    return df.apply( lambda r: peptides.getAAFraction(str(r.peptide), amino_acids), 1)

def netCharge(df):
    """Net peptide charge for dataframe with sequences"""
    return df.apply( lambda r: peptides.netCharge(r.peptide),1)

def isoelectricPoint(df):
    def getpi(seq):
        X = ProteinAnalysis(seq)
        return X.isoelectric_point()
    return df.apply( lambda r: getpi(r.peptide),1)

def _center_nmer(x, n):
    """Get n-mer sequence for a peptide centered in the middle.
    This should be applied to a dataframe per row."""

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

def _split_nmer(x, n, key):
    """Row based method to split a peptide in to multiple nmers
    if it's too large"""

    size = x.end-x.start
    m = 3 #margin
    if size <= n+m:
        return pd.Series(_center_nmer(x, n))
    else:
        o=size%n
        #print (x.peptide, size, o)
        seqs=[]
        l = int(math.ceil(float(size)/n))
        if x.start==0: s=1
        else: s=0
        for i in range(s, size, n):
            if i+n>size:
                seqs.append(x[key][o:o+n])
            else:
                seqs.append(x[key][i:i+n])
        seqs = pd.Series(seqs)
        return seqs

def getNmer(df, genome, length=20, seqkey='peptide', how='center'):
    """
    Get n-mer peptide surrounding a set of sequences using the host
    protein sequence.
    Args:
        df: input dataframe with sequences
        genome: genome dataframe with host sequences
        length: length of nmer to return
        seqkey: column name of sequence to be processed
        how: method to create the n-mer, split will try to split up
            the sequence into overlapping n-mes of length is larger than size
    Returns:
        pandas Series with nmer values
    """

    temp = df.merge(genome[['locus_tag','gene','translation','length']],
                    left_on='name',right_on='locus_tag',how='left')

    if not 'end' in list(temp.columns):
        temp = base.getCoords(temp)
    temp =  base.getCoords(temp)
    if how == 'center':
        res = temp.apply( lambda r: _center_nmer(r, length), 1)
    elif how == 'split':
        res = temp.apply( lambda r: _split_nmer(r, length, seqkey), 1)
        res.index = temp.index
        res = res.stack().drop_duplicates()
        res.index = res.index.droplevel(1)
    return res

def getOverlaps(binders1, binders2, label='overlaps', how='inside'):
    """
    Overlaps for 2 sets of sequences where positions in host sequence are stored
    in each dataframe.
    Args:
        binders1: first set of sequences, a dataframe with pos field
        binders2: second set of sequences
        label: label for overlaps column
        how: may be 'any' or 'inside'
    Returns:
        First DataFrame with no. of overlaps stored in a new column
    """
    new=[]
    a = base.getCoords(binders1)
    b = base.getCoords(binders2)

    def overlap(x,y):
        f=0
        #print x['name'],x.peptide
        #print x.start,x.end
        for i,r in y.iterrows():
            if how == 'inside':
                if ((x.start<=r.start) & (x.end>=r.end)):
                    f+=1
            elif how == 'any':
                if ((x.start<r.start) & (x.end>r.start)) or \
                   ((x.start>r.start) & (x.start<r.end)):
                    f+=1
            #print r.start,r.end, f
        return f

    for n,df in a.groupby('name'):
        found = b[b.name==n]
        df[label] = df.apply(lambda r: overlap(r,found),axis=1)
        new.append(df)
    result = pd.concat(new)
    #print (len(a), len(b))
    print ('%s with overlapping binders' %len(result[result[label]>0]))
    return result

def getLocalOrthologs(seq, db):
    """Get alignment for a protein using local blast db"""


    result_handle = open("my_blast.xml")
    df = sequtils.getBlastResults(result_handle)
    return df

def getOrthologs(seq, db=None, expect=10, hitlist_size=400, equery=None,
                 email=''):
    """
    Fetch orthologous sequences using online blast and return the records
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
        sequtils.localBlast(db, 'tempseq.faa', output='my_blast.xml', maxseqs=30)
        result_handle = open("my_blast.xml")
        df = sequtils.getBlastResults(result_handle)
    else:
        try:
            result_handle = NCBIWWW.qblast("blastp", "nr", seq, expect=expect,
                                  hitlist_size=500,entrez_query=equery)
            time.sleep(2)
            savefile = open("my_blast.xml", "w")
            savefile.write(result_handle.read())
            savefile.close()
            result_handle = open("my_blast.xml")
            df = sequtils.getBlastResults(result_handle, local=False)
        except Exception as e:
            print ('blast timeout')
            return

    df = df.drop(['subj','positive','query_length','score'],1)
    df.drop_duplicates(subset=['definition','perc_ident'], inplace=True)
    df = df[df['perc_ident']!=100]
    return df

def alignBlastResults(df, aln=None, idkey='accession', productkey='definition'):
    """
    Get gapped alignment from blast results using muscle aligner.
    """

    sequtils.dataframe2Fasta(df, idkey=idkey, seqkey='sequence',
                        descrkey=productkey, outfile='blast_found.faa')
    aln = sequtils.muscleAlignment("blast_found.faa")
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
    try:
        ind = alnrows.species #alnrows.accession
    except:
        ind = alnrows.name
    s = pd.DataFrame(f,columns=seqs,index=ind)
    s = s.replace(-1,np.nan)
    s[s>0] = 1
    res = s.count()
    return s

def epitopeConservation(peptides, alnrows=None, proteinseq=None, blastresult=None,
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
        Matrix of
    """

    import seaborn as sns
    sns.set_context("notebook", font_scale=1.4)

    if alnrows is None:
        if proteinseq == None:
            print ('protein sequence to blast or alignment required')
            return
        if blastresult == None or not os.path.exists(blastresult):
            blr = getOrthologs(proteinseq, equery=equery, blastdb=blastdb)
            if blr is None:
                return
            #if filename == None: filename = 'blast_%s.csv' %label
            blr.to_csv(blastresult)
        else:
            blr = pd.read_csv(blastresult, index_col=0)
        #blr = blr[blr.perc_ident>=perc_ident]
        alnrows, aln = alignBlastResults(blr)
        #print (sequtils.formatAlignment(aln))

    if 'perc_ident' in alnrows.columns:
        alnrows = alnrows[alnrows.perc_ident>=perc_ident]
    if 'definition' in alnrows.columns:
        alnrows['species'] = alnrows.definition.apply(get_species_name)
    c = find_conserved_sequences(peptides, alnrows).T
    #print (c)
    c = c.dropna(how='all')
    c = c.reindex_axis(c.sum(1).sort_values().index)
    if len(c) == 0:
        print ('no conserved epitopes in any sequence')
        return
    #c['fraction'] = (c.count(1)/len(c.columns)).round(2)
    h=int(len(c)/2.)+2
    sns.heatmap(c, linewidths=1, cmap='summer', linecolor='black',
                cbar=False, square=True)
    return c

'''def setBlastLink(df):
    def makelink(x):
        return '<a href=http://www.ncbi.nlm.nih.gov/protein/%s> %s </a>' %(x,x)
    df['accession'] = df.accession.apply(makelink)
    return df'''

def findClusters(binders, dist=None, min_binders=2, min_size=12, max_size=50,
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
    length = len(binders.head(1).peptide.max())
    if dist == None:
        dist = length+1
        print ('using dist for clusters: %s' %dist)
    for n,b in grps:
        if len(b)==0: continue
        clusts = base.dbscan(b,dist=dist,minsize=min_binders)
        if len(clusts) == 0:
            continue
        for c in clusts:
            gaps = [c[i]-c[i-1] for i in range(1,len(c))]
            C.append([n,min(c),max(c)+length,len(c)])

    if len(C)==0:
        print ('no clusters')
        return pd.DataFrame()
    x = pd.DataFrame(C,columns=['name','start','end','binders'])
    x['clustersize'] = (x.end-x.start)
    x = x[x.clustersize>=min_size]
    x = x[x.clustersize<=max_size]

    #if genome data available merge to get peptide seq
    if genome is not None:
        x = x.merge(genome[['locus_tag','translation','gene']],
                    left_on='name',right_on='locus_tag')
        x[colname] = x.apply(lambda r: r.translation[r.start:r.end], 1)
        x = x.drop(['locus_tag','translation'],1)
    x = x.sort_values(by=['binders'],ascending=False)
    x = x.reset_index(drop=True)
    print ('%s clusters found in %s proteins' %(len(x),len(x.groupby('name'))))
    print
    return x

def randomizedLists(df, n=94, seed=8, filename='peptide_lists'):
    """
    Return a randomized lists of sequences from a dataframe. Used for
    providing peptide lists for aassying etc.
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

def genomeAnalysis(datadir,label,gname,method):
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
    return res

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

def getSeqDepot(seq):
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

def testFeatures():
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
    #b = P.getBinders(data=df)
    #print b[:20]
    base.getScoreDistributions(method, path)
    return

def testBcell(gname):
    path='test'
    gfile = os.path.join(genomespath,'%s.gb' %gname)
    df = sequtils.genbank2Dataframe(gfile, cds=True)
    names=['VP24']
    P = base.getPredictor('bcell')
    P.iedbmethod='Chou-Fasman'
    P.predictProteins(df,names=names,save=True,path=path)
    print (P.data)
    return

def testconservation(label,gname):
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
    b = P.getBinders()'''
    return

def findConservedPeptide(peptide, recs):
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

def getPredictions(path,tag,method='tepitope',q=0.96):
    """Get predictions from file system"""

    q=round(q,2)
    #preds = OrderedDict()
    cutoffs = {}
    filename = os.path.join(path, tag+'.mpk')
    if not os.path.exists(filename):
        return
    df = pd.read_msgpack(filename)
    pred = base.getPredictor(name=method, data=df)
    cutoffs = pred.allelecutoffs = getCutoffs(path, method, q)
    pred = pred
    return pred

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
