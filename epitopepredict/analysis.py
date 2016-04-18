#!/usr/bin/env python

"""
    epitopepredict analysis methods for workflows
    Created September 2013
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, shutil, string, types
import csv, glob, pickle, itertools
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
from . import base, sequtils, tepitope, utilities

home = os.path.expanduser("~")
#fix paths!
genomespath = os.path.join(home, 'epitopedata')
datadir = os.path.join(home, 'testpredictions')

def plotheatmap(df, ax=None, cmap='Blues'):

    if ax==None:
        fig=plt.figure()
        ax=fig.add_subplot(111)
    else:
        fig = ax.get_figure()
    df = df._get_numeric_data()
    hm=ax.pcolor(df,cmap=cmap)
    #fig.colorbar(hm, ax=ax)
    ax.set_xticks(np.arange(0.5, len(df.columns)))
    ax.set_yticks(np.arange(0.5, len(df.index)))
    ax.set_xticklabels(df.columns, minor=False, fontsize=10,rotation=45)
    ax.set_yticklabels(df.index, minor=False, fontsize=8)
    ax.set_ylim(0, len(df.index))
    hm.set_clim(0,1)
    plt.tight_layout()
    return

'''def compareBindingData(exp, pred, seqkey, datakey, allele):
    """Compare to experimental binding data in a csv file.
     pred: the predictor
     seqkey: col name for sequences
     datakey col name for exp binding/score value"""

    peptides = list(exp[seqkey])
    pred.predict(peptides=peptides, allele=allele)
    scorekey = pred.scorekey
    data = pd.read_csv(os.path.join(tepitope.tepitopedir,
                        'predictions/HLA-DRB1*0101.csv'))
    pred.cutoff=0
    merged = pd.merge(pred.data, exp,
                  left_on='peptide', right_on=seqkey,
                  suffixes=['_1', '_2'])
    def cutoff(x,c,func):
        if func(x,c): return 1
        else: return 0

    merged['pred'] = merged[scorekey].apply(lambda x:cutoff(x,pred.cutoff,op.gt))
    merged['exp'] = merged[datakey].apply(lambda x:cutoff(x,.426,op.gt))
    #print merged['exp']
    P=merged[scorekey]
    E=merged[datakey]
    from sklearn import metrics
    auc = metrics.roc_auc_score(merged['exp'],merged['pred'])

    print '%s, %s samples' %(allele,len(E))
    print auc
    print Predict.aucCrossValidation(merged['exp'].values,merged['pred'].values)
    return

def compareProteins(df, pred, exp, scorefield):
    """Compare whole antigen scores to % best binders per protein"""

    binders=[]
    names = []
    cds = df[df.type=='CDS']
    for i,row in list(cds.iterrows())[:10]:
        seq = row['translation']
        name = row['locus_tag']
        print name,seq
        names.append(name)
        pred.predict(sequence=seq)
        binders.append(pred.getBinders(method='cutoff'))
    return

def getMatchingPredictions(pred1, pred2, method='cutoff'):
    """Compare 2 predictors binders"""

    data1 = pred1.getBinders(method=method)
    data2 = pred2.getBinders(method=method)
    #data1 = pred1.data
    #data2 = pred2.data
    merged = pd.merge(data1, data2,
                      left_on=['peptide','name'], right_on=['peptide','name'],
                      suffixes=['_1', '_2'])
    union = pd.merge(data1, data2, how='outer',
                      left_on=['peptide','name'], right_on=['peptide','name'])

    #merged = merged[merged.core_1==merged.core_2]
    #print merged[:20]
    k1 = pred1.scorekey; k2 = pred2.scorekey
    if k1 == k2:
        k1 += '_1'
        k2 += '_2'
    x=merged[k1]
    y=merged[k2]

    print '%s/%s shared binders' %(len(merged),len(union))
    return merged, x, y

def comparePredictors(pred1, pred2,
                      names=None, allele='HLA-DRB1*0101'):
    """Compare 2 predictors with various metrics and plot output.
       Input is a dataframe with sequence records"""

    from matplotlib_venn import venn2
    import pylab as plt
    f = plt.figure(figsize=(10,10))
    ax = f.add_subplot(221)

    binders1 = pred1.getBinders('cutoff')
    binders2 = pred2.getBinders('cutoff')
    names = dict(list(pred1.data.groupby('name'))).keys()

    #get merged binders
    m,x,y = getMatchingPredictions(pred1, pred2, method='cutoff')

    ax.plot(x,y,'o',ms=3,alpha=0.5)
    ax.set_xlabel(pred1.name)
    ax.set_ylabel(pred2.name)
    ind=np.arange(len(names))
    b1 = list(binders1.peptide)
    b2 = list(binders2.peptide)
    #print list(set(names1) & set(names2))
    groups1 = dict(list(binders1.groupby('name')))
    groups2 = dict(list(binders2.groupby('name')))
    prots1 = groups1.keys()
    prots2 = groups2.keys()

    ax1 = f.add_subplot(222)
    ax1.set_title('proteins overlap')
    venn2([set(prots1), set(prots2)], set_labels = (pred1.name,pred2.name))
    ax3=f.add_subplot(212)
    venn2([set(b1), set(b2)], set_labels = (pred1.name,pred2.name))
    f.suptitle('%s vs. %s' %(pred1.name,pred2.name),fontsize=16)
    plt.show()
    return'''

def getNmer(df, n=20, key='translation'):
    """Get 20mer peptide"""

    def getseq(x):
        n=20
        size=len(x[key])
        if size<n:
            o = int((n-size)/2.0)+1
            s = x[key][x.start-o:x.end+o][:20]
        else:
            s = x[key][x.start:x.end]
        return s
    x = df.apply(getseq,1)
    return x

def getOverlappingBinders(binders1, binders2, label='overlap'):
    """Overlap for binders with any set of peptides with start/end cols"""

    new=[]
    def overlap(x,b):
        f = b[(b.pos>x.start) & (b.pos<x.end)]
        #print x.locus_tag,x.start,x.end,x.peptide,len(f) #,f.peptide,f.pos
        return len(f)
    for n,df in binders1.groupby('name'):
        b = binders2[binders2.name==n]
        df[label] = df.apply(lambda r: overlap(r,b),axis=1)
        new.append(df)
    result = pd.concat(new)
    print ('%s with overlapping binders' %len(result[result[label]>0]))
    return result

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
        print ('running blast..')
        result_handle = NCBIWWW.qblast("blastp", "nr", seq, expect=expect,
                              hitlist_size=500,entrez_query=equery)
        time.sleep(2)
    except:
        print ('blast timeout')
        return
    savefile = open("my_blast.xml", "w")
    savefile.write(result_handle.read())
    savefile.close()
    result_handle = open("my_blast.xml")

    df = sequtils.getBlastResults(result_handle)
    df['accession'] = df.subj.apply(lambda x: x.split('|')[3])
    df['definition'] = df.subj.apply(lambda x: x.split('|')[4])
    df = df.drop(['subj','positive','query_length','score'],1)
    print (len(df))
    df.drop_duplicates(subset=['definition'], inplace=True)
    df = df[df['perc_ident']!=100]
    print (len(df))
    #df = getAlignedBlastResults(df)
    return df

def getAlignedBlastResults(df,aln=None,idkey='accession',productkey='definition'):
    """Get gapped alignment from blast results """

    sequtils.dataframe2Fasta(df, idkey=idkey, seqkey='sequence',
                        productkey=productkey, outfile='blast_found.faa')
    aln = sequtils.muscleAlignment("blast_found.faa")
    alnrows = [[a.id,str(a.seq)] for a in aln]
    alndf = pd.DataFrame(alnrows,columns=['accession','seq'])
    #res = df.merge(alndf, left_index=True, right_index=True)
    res = df.merge(alndf, on=['accession'])
    res = res.drop('sequence',1)
    #get rid of duplicate hits
    #res.drop_duplicates(subset=['definition','seq'], inplace=True)
    res = res.sort('identity',ascending=False)
    print ('%s hits, %s filtered' %(len(df), len(res)))
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

def findClusters(binders, method, dist=None, minsize=3,
                 genome=None):
    """Get clusters of binders for all predictions"""

    C=[]
    grps = list(binders.groupby('name'))
    print ('%s proteins with binders' %len(grps))
    length = len(binders.head(1).peptide.max())
    if dist==None:
        dist = length+1
        print ('using dist for clusters: %s' %dist)
    for n,b in grps:
        if len(b)==0: continue
        clusts = base.dbscan(b,dist=dist,minsize=minsize)
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
    x['density'] = x.binders/(x.end-x.start)
    x['method'] = method

    if genome is not None:
        temp = x.merge(genome[['locus_tag','gene','translation']],
                    left_on='name',right_on='locus_tag')
        x['peptide'] = getNmer(temp)
    x = x.sort_values(by=['binders','density'],ascending=False)
    print ('%s clusters found in %s proteins' %(len(x),len(x.groupby('name'))))
    print
    return x

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

def getLocalOrthologs(seq, db):
    """Get alignment for a protein using local blast db"""

    SeqIO.write(SeqRecord(Seq(seq)), 'tempseq.faa', "fasta")
    sequtils.doLocalBlast(db, 'tempseq.faa', output='my_blast.xml', maxseqs=30)
    result_handle = open("my_blast.xml")
    df = sequtils.getBlastResults(result_handle)
    return df

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
