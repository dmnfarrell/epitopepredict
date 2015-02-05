#!/usr/bin/env python

"""
    MHC Epitope analysis
    Created September 2013
    Copyright (C) Damien Farrell
"""

import sys, os, shutil, string, types
import csv, glob, pickle, itertools
import re
import time, random
import collections
from operator import itemgetter
#import matplotlib
#matplotlib.use('agg')
import pylab as plt
import numpy as np
import pandas as pd
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import base, sequtils, tepitope, utilities

home = os.path.expanduser("~")
genomespath = os.path.join(home, 'epitopedata')
datadir = os.path.join(home, 'testpredictions')

def getAllBinders(path, method='tepitope', n=3, cutoff=0.95):
    """Get all binders from a set of proteins in path"""

    print 'getting binders..'
    binders = []
    m=method
    if m=='bcell': return #not applicable
    l=9
    P = base.getPredictor(m)
    files = glob.glob(os.path.join(path, '*.mpk'))
    #get allele specific cutoffs
    P.allelecutoffs = getCutoffs(path, method, cutoff, overwrite=True)
    for f in files:
        df = pd.read_msgpack(f)
        b = P.getPromiscuousBinders(data=df,n=n)
        #print b[:5]
        binders.append(b)
    result = pd.concat(binders)
    return result

def getCutoffs(path, method, q=0.98, overwrite=False):
    """Get global cutoffs for predictions in path"""

    quantfile = os.path.join(path,'quantiles.csv')
    if not os.path.exists(quantfile) or overwrite==True:
        base.getScoreDistributions(method, path)
    quantiles = pd.read_csv(quantfile,index_col=0)
    cutoffs = dict(quantiles.ix[q])
    return cutoffs

def getOrthologs(seq,expect=10,hitlist_size=400,equery=None):
    """Fetch orthologous sequences using blast and return the records
        as a dataframe"""

    from Bio.Blast import NCBIXML,NCBIWWW
    from Bio import Entrez, SeqIO
    Entrez.email = "anon.user@ucd.ie"
    #entrez_query = "mycobacterium[orgn]"
    #db = '/local/blast/nr'
    #SeqIO.write(SeqRecord(Seq(seq)), 'tempseq.faa', "fasta")
    #sequtils.doLocalBlast(db, 'tempseq.faa', output='my_blast.xml', maxseqs=100, evalue=expect)

    try:
        print 'running blast..'
        result_handle = NCBIWWW.qblast("blastp", "nr", seq, expect=expect,
                              hitlist_size=500,entrez_query=equery)
        time.sleep(2)
    except:
        print 'blast timeout'
        return
    savefile = open("my_blast.xml", "w")
    savefile.write(result_handle.read())
    savefile.close()
    result_handle = open("my_blast.xml")

    df = sequtils.getBlastResults(result_handle)
    df['accession'] = df.subj.apply(lambda x: x.split('|')[3])
    df['definition'] = df.subj.apply(lambda x: x.split('|')[4])
    df = df.drop(['subj','positive','query_length','score'],1)
    print len(df)
    df.drop_duplicates(subset=['definition'], inplace=True)
    df = df[df['perc_ident']!=100]
    print len(df)
    #df = getAlignedBlastResults(df)
    return df

def getAlignedBlastResults(df,aln=None):
    """Get gapped alignment from blast records"""

    sequtils.dataframe2Fasta(df, idkey='accession', seqkey='sequence',
                        productkey='definition', outfile='blast_found.faa')
    aln = sequtils.muscleAlignment("blast_found.faa")
    alnrows = [[a.id,str(a.seq)] for a in aln]
    alndf = pd.DataFrame(alnrows,columns=['accession','seq'])
    res = df.merge(alndf, on='accession')
    res = res.drop('sequence',1)
    #get rid of duplicate hits
    #res.drop_duplicates(subset=['definition','seq'], inplace=True)
    res = res.sort('identity',ascending=False)
    print '%s hits, %s filtered' %(len(df), len(res))
    print res[['accession','definition']][-20:]
    return res

def setBlastLink(df):
    def makelink(x):
        return '<a href=http://www.ncbi.nlm.nih.gov/protein/%s> %s </a>' %(x,x)
    df['accession'] = df.accession.apply(makelink)
    return df

def alignment2Dataframe(aln):
    """Blast results alignment 2 dataframe for making tables"""

    alnrows = [[a.id,str(a.seq)] for a in aln]
    df = pd.DataFrame(alnrows,columns=['name','seq'])
    return df

def findClusters(binders, method, dist=None, minsize=3):
    """Get clusters of binders for all predictions"""

    C=[]
    grps = list(binders.groupby('name'))
    print '%s proteins with binders' %len(grps)
    length = len(binders.head(1).peptide.max())
    if dist==None:
        dist = length+1
        print 'using dist for clusters: %s' %dist
    for n,b in grps:
        if len(b)==0: continue
        clusts = base.dbscan(b,dist=dist,minsize=minsize)
        if len(clusts) == 0:
            continue
        for c in clusts:
            gaps = [c[i]-c[i-1] for i in range(1,len(c))]
            C.append([n,min(c),max(c)+length,len(c)])

    if len(C)==0:
        print 'no clusters'
        return
    x = pd.DataFrame(C,columns=['name','start','end','binders'])
    x['clustersize'] = (x.end-x.start)
    x['density'] = x.binders/(x.end-x.start)
    x['method'] = method
    x = x.sort(['binders','density'],ascending=False)
    print '%s clusters found in %s proteins' %(len(x),len(x.groupby('name')))
    print
    return x

def testgenomeanalysis(label,gname,method):
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

    #fig=plt.figure(figsize=(6,6))
    #ax=fig.add_subplot(111)
    #b.hist(P.scorekey,bins=30,alpha=0.8,ax=ax)
    #res.sort('order').plot(kind='bar',x='order',y='perc',ax=ax)
    #plt.show()
    top = b.groupby('peptide').agg({P.scorekey:np.mean,'allele':np.max,
                    'name': lambda x: x}).reset_index()
    top = top.sort(P.scorekey,ascending=P.rankascending)
    print b
    cl = findClusters(b, method, dist=9, minsize=3)
    if cl is not None:
        gc = cl.groupby('name').agg({'density':np.max})
        res = res.merge(gc,left_on='locus_tag',right_index=True)
    print res[:10]

    return

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
    print fastafmt
    print row.to_dict()
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
                        save=True,path=path)
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
    print P.data
    return

def testconservation(label,gname):
    """Conservation analysis"""

    tag='VP24'
    pd.set_option('max_colwidth', 800)
    gfile = os.path.join(genomespath,'%s.gb' %gname)
    g = sequtils.genbank2Dataframe(gfile, cds=True)
    res = g[g['locus_tag']==tag]
    seq = res.translation.head(1).squeeze()
    print seq
    #alnrows = getOrthologs(seq)
    #alnrows.to_csv('blast_%s.csv' %tag)
    alnrows = pd.read_csv('blast_%s.csv' %tag,index_col=0)
    alnrows.drop_duplicates(subset=['accession'], inplace=True)
    alnrows = alnrows[alnrows['perc_ident']>=60]
    seqs=[SeqRecord(Seq(a.sequence),a.accession) for i,a in alnrows.iterrows()]
    print seqs[:2]
    sequtils.distanceTree(seqs=seqs)#,ref=seqs[0])
    #sequtils.ETETree(seqs, ref, metric)
    #df = sequtils.getFastaProteins("blast_found.faa",idindex=3)
    '''method='tepitope'
    P = base.getPredictor(method)
    P.predictSequences(df,seqkey='sequence')
    b = P.getBinders()'''
    return

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

