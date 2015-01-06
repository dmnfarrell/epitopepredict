#!/usr/bin/env python

"""
    MHC prediction base module for core classes
    Created November 2013
    Copyright (C) Damien Farrell
"""

import sys, os, shutil, string
import csv, glob, pickle
import time, StringIO
import operator as op
import re, types
import subprocess
import numpy as np
import pandas as pd
import pylab as plt
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from matplotlib.ticker import MaxNLocator
import utilities, peptides, genome, threading, tepitope

home = os.path.expanduser("~")
datadir = os.path.join(home, 'mhcdata')
iedbmethods = ['arbpython','comblib','consensus3','IEDB_recommended',
               'NetMHCIIpan','nn_align','smm_align','tepitope']
iedbsettings = {'cutoff_type': 'none', 'pred_method': 'IEDB_recommended',
            'output_format': 'ascii', 'sort_output': 'position_in_sequence',
            'sequence_format': 'auto', 'allele': 'HLA-DRB1*01:01', 'length':'11',
            'sequence_file': None}
iedbkeys = {'consensus3': ['Allele','Start','End','Sequence','consensus_percentile',
            'comblib_core','comblib_score','comblib_percentile','smm_core','smm_score',
            'smm_percentile','nn_core','nn_score','nn_percentile','Sturniolo core',
            'Sturniolo score','Sturniolo percentile'],
        'IEDB_recommended': ['Allele','Start','End','Sequence','consensus_percentile',
            'comblib_core','comblib_score','comblib_percentile','smm_core','smm_score',
            'smm_percentile','nn_core','nn_score','nn_percentile','netMHCIIpan_core',
            'netMHCIIpan_score','netMHCIIpan_percentile','Sturniolo core',
            'Sturniolo score','Sturniolo percentile','methods'],
        'NetMHCIIpan': ['Allele','Start','End','Core','Sequence','IC50']}

#based on MTB proteome score dists, need to make this more general..
globalcutoffs = {'tepitope': {9: {'HLA-DRB1*0801': 1.8, 'HLA-DRB3*0201': 1.72,
                'HLA-DRB1*1101': 1.2, 'HLA-DRB1*1401': 1.83, 'HLA-DRB1*0301': 2.96,
                'HLA-DRB1*0401': 1.4, 'HLA-DRB3*0101': 1.82, 'HLA-DRB1*1301': 2.21,
                'BoLA-DRB3*0401': 2.2, 'BoLA-DRB3*2101': 1.56, 'BoLA-DRB3*3701': 1.63,
                'BoLA-DRB3*4101': 2.0, 'BoLA-DRB3*6301': 2.09, 'BoLA-DRB3*2002': 1.67,
                'BoLA-DRB3*3001': 1.95, 'BoLA-DRB3*1601': 1.27, 'BoLA-DRB3*1901': 2.25,
                'BoLA-DRB3*1201': 2.12}},
                 'netmhciipan': {9: {'DRB1*0801': 0.29, 'DRB1*1301': 0.16, 'DRB3*0101': 0.09,
                  'DRB3*0201': 0.19, 'DRB1*0301': 0.08, 'DRB1*1101': 0.24,
                  'DRB1*1401': 0.13, 'DRB1*0401': 0.11},
                  11: {'DRB1*0801': 0.54, 'DRB1*1301': 0.41, 'DRB3*0101': 0.27,
                 'DRB3*0201': 0.49, 'DRB1*0301': 0.29, 'DRB1*1101': 0.5,
                 'DRB1*1401': 0.32, 'DRB1*0401': 0.3},
                  15: {'DRB1*0801': 0.7, 'DRB1*1301': 0.66, 'DRB3*0101': 0.53,
                 'DRB3*0201': 0.75, 'DRB1*0301': 0.65, 'DRB1*1101': 0.67,
                 'DRB1*1401': 0.51, 'DRB1*0401': 0.55}}
                 }

def first(x):
    return x.iloc[0]

def getIEDBRequest(seq, alleles='HLA-DRB1*01:01', method='consensus3'):
    import requests
    url = 'http://tools.iedb.org/tools_api/mhcii/'
    values = {'method' : method,
              'sequence_text' : seq,
              'allele' : alleles }
    r=requests.post(url,data=values)
    df=pd.read_csv(StringIO.StringIO(r.content),sep='\t')
    #df=df.drop(['nn_align_core','nn_align_ic50','nn_align_rank'])
    return df

def venndiagram(names,labels,ax=None):
    from matplotlib_venn import venn2,venn3
    f=None
    if ax==None:
        f=plt.figure(figsize=(4,4))
        ax=f.add_subplot(111)
    if len(names)==2:
        n1,n2=names
        v = venn2([set(n1), set(n2)], set_labels=labels)
    elif len(names)==3:
        n1,n2,n3=names
        v = venn3([set(n1), set(n2), set(n3)], set_labels=labels)
    ax.axis('off')
    #f.patch.set_visible(False)
    ax.set_axis_off()
    #plt.tight_layout()
    return f

def getOverlapping(index, s, length=9, cutoff=25):
    """Get all mutually overlapping kmers within a cutoff area"""
    g=[s]
    vals = [i for i in range(s, s+cutoff) if i in index]
    for i in range(len(vals)-1):
        if vals[i+1]<=vals[i]+length:
            g.append(vals[i+1])
        else:
            break
    return g

def checkMembers(g,clusts):
    """Check if a group intersects any of the current clusters"""
    for i in clusts:
        common = list(set(g) & set(i))
        if len(common)>0 or len(g)<2:
            #print i,common
            return False
    return True

def getClusters(B, clustlen=25, cutoff=0.05):
    """Get clusters of binders from a set of predictions
      df: a pandas dataframe with one set of predictions per row"""

    nmer = len(B.iloc[0].peptide)
    overlap = clustlen - nmer
    #print clustlen, nmer, overlap
    locs = pd.Series(B.peptide.values,index=B.pos).to_dict()
    #ad hoc method to get overlapping epitopes and
    #return largest unique groups as clusters
    groups=[]
    for i in locs:
        g = getOverlapping(locs, int(i), overlap, clustlen)
        groups.append(g)
    ranked = sorted(groups, key=len, reverse=True)
    clusts=[]
    for g in ranked:
        if checkMembers(g, clusts) == True:
            clusts.append(g)
    #print clusts
    return clusts

def dbscan(B=None, x=None, dist=7, minsize=4):
    """Density-Based Spatial clustering. Finds core samples of
      high density and expands clusters from them."""
    from sklearn.cluster import DBSCAN
    if B is not None:
        if len(B)==0:
            return
        x = sorted(B.pos.astype('int'))
    X = np.array(zip(x,np.zeros(len(x))), dtype=np.int)
    db = DBSCAN(eps=dist, min_samples=minsize)
    db.fit(X)
    labels = db.labels_
    n_clusters_ = len(set(labels))
    clusts=[]
    for k in range(n_clusters_):
        my_members = labels == k
        #print "cluster {0}: {1}".format(k, X[my_members, 0])
        if len(X[my_members, 0])>0:
            clusts.append(list(X[my_members, 0]))
    #print clusts
    return clusts

def getPredictor(name='tepitope', **kwargs):
    if name == 'netmhciipan':
        return NetMHCIIPanPredictor(**kwargs)
    elif name == 'iedbmhc1':
        return IEDBMHCIPredictor(**kwargs)
    elif name == 'iedbmhc2':
        return IEDBMHCIIPredictor(**kwargs)
    elif name == 'tepitope':
        return TEpitopePredictor(**kwargs)
    elif name == 'threading':
        return ThreadingPredictor(**kwargs)
    elif name == 'modeller':
        return ModellerPredictor(**kwargs)
    else:
        print 'no such predictor %s' %name
        return

def createTempSeqfile(sequences, seqfile='tempseq.fa'):
    if type(sequences) is types.StringType:
        sequences=[sequences]
    out = open(seqfile, 'w')
    i=1
    for seq in sequences:
        SeqIO.write(SeqRecord(Seq(seq),id='temp%s'%i,
                    description='temp'), out, 'fasta')
        i+=1
    out.close()
    return seqfile

def getSequence(seqfile):
    """Get sequence from fasta file"""
    recs = list(SeqIO.parse(seqfile, 'fasta'))[0]
    sequence = recs.seq.tostring()
    return sequence

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

def tmhmm(infile=None,genome=None):
    """Get TMhmm predictions"""
    if genome != None:
        seqfile = Genome.genbank2Fasta(genome)
        tempfile = 'tmhmm_temp.txt'
        cmd = 'tmhmm %s > %s' %(seqfile,tempfile)
        infile = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    tmpred = pd.read_csv(infile,delim_whitespace=True,comment='#',
                      names=['locus_tag','v','status','start','end'])
    tmpred = tmpred.dropna()
    print 'tmhmm predictions for %s proteins' %len(tmpred.groupby('locus_tag'))
    lengths=[]
    for i,row in tmpred.iterrows():
        if row.status == 'TMhelix':
            lengths.append(row.end-row.start)
    #print np.mean(lengths), np.std(lengths)
    return tmpred

def netSurfP(infile=None, genome=None):
    """NetsurfP predictions"""
    if genome != None:
        seqfile = Genome.genbank2Fasta(infile)
        tempfile = 'netsurfp_temp.txt'
        cmd = 'netsurfp -i %s > %s' %(seqfile,tempfile)
        infile = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    names = ['class','aa','locus_tag', 'pos','rsa', 'asa','zscore']
    pred = pd.read_csv(infile,delim_whitespace=True,comment='#',
                        names=names)
    pred = pred.dropna()
    return pred

def getTMRegions(tm, name):
    df = tm[tm.locus_tag==name]
    df = df[df.status=='TMhelix']
    return zip(df.start,df.end)

def getOverlappingBinders(B):
    for df in B:
        df.set_index('peptide',inplace=True)
    print B
    x = pd.concat(B,join='inner')
    #x =  x[['pos']].sort('pos')
    return x

def compareBindingData(exp, pred, seqkey, datakey, allele):
    """Compare to experimental binding data in a csv file.
     pred: the predictor
     seqkey: col name for sequences
     datakey col name for exp binding/score value"""

    peptides = list(exp[seqkey])
    Tepitope.alpha=10
    pred.predict(peptides=peptides, allele=allele)
    scorekey = pred.scorekey
    data = pd.read_csv(os.path.join(Tepitope.tepitopedir,
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
    return

def plotBars(x, y, ind, ax, label=None,
                color='b', ylim=None, xlabels=False):
    """Bar plot showing scores over an index"""

    bars = ax.bar(x,y,alpha=0.8,color=color,align='center',lw=0)
    ax.set_xlim([min(ind) - 0.5, max(ind) + 0.5])
    if label!=None:
        ax.text(0.8,0.6,label,transform=ax.transAxes,fontsize=14)
    if ylim!=None:
        ax.set_ylim(ylim)
    if xlabels == True:
        ax.set_xticks(ind+0.5)
        ax.set_xticklabels(df.index.tolist(),rotation=45)
    return bars

def plotGroups(grps, ycol, sortcol='pos', ylim=None,
                xlabels=False, kind='bar', ax=None, label=None):
    """Plot multiple dataframe groups"""

    i=0
    size=len(grps)
    cmap = plt.cm.get_cmap('jet')
    fig=None
    if ax == None:
        if size == 1:
            fig = plt.figure()
            grid = [fig.add_subplot(111)]
        else:
            fig, grid = plt.subplots(nrows=size, ncols=1, figsize=(12,6),
                                     sharex=True, sharey=True)
    currax=ax
    for a, df in grps:
        if ax==None:
            currax = grid[i]
        c = cmap(float(i)/(size))
        if label == None: l=a
        else: l=label
        plotScores(df,ycol,currax,sortcol=sortcol,label=l,color=c,
                   ylim=ylim,xlabels=xlabels,kind=kind)
        currax.yaxis.set_major_locator(MaxNLocator(3))
        i+=1
    return fig

def plotPerAllele(P, names=None, path='plots', kind='bar'):
    """Plot predictions for each protein per allele"""

    df = P.data
    scorekey = P.scorekey
    for name,protdf in df.groupby('name'):
        if names != None and name not in names:
            continue
        grps = protdf.groupby('allele')
        if P.cutoff < 0:
            highest = min(protdf[scorekey])
        else:
            highest = max(protdf[scorekey])
        lims = (P.cutoff,highest)
        f = plotGroups(grps,scorekey,ylim=lims,kind=kind)
        plt.tight_layout(h_pad=0.5)
        #f.suptitle(name,fontsize=16)
        filename = os.path.join(path,name+'_%s.png' %P.name)
        f.savefig(filename, dpi=80)
    return filename

def plotRegions(locs, ax,color='g',**kwargs):
    """Draw boxes to highlight areas on a plotted sequence"""

    if locs==None:
        return
    from matplotlib.patches import Rectangle
    for i in locs:
        x1=min(i)-0.5; x2=max(i)-0.5
        y1,y2 = ax.get_ylim()
        ax.add_patch(Rectangle((x1, y1), x2-x1+1, y2,color=color,
                    **kwargs))
    return

def plotHeatMap(pred, name):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    p = pred.reshape(name)
    p=p.sort('pos')
    x = p.drop(['peptide','pos','mean'],1)
    r = x.transpose()
    plt.pcolor(r)
    return

def plotMultipleScores(predictors, names, **kwargs):
    """Plot multiple proteins"""

    for n in names:
        print n
        f=plotPredictorScores(predictors, n, **kwargs)
    return

def plotPredictorScores(predictors, name, cldist=7, n=2,
                        path='plots', save=True, dpi=80, **kwargs):
    """Plot multiple predictors in same plot"""

    size = len(predictors)
    if size == 1:
        fig = plt.figure(figsize=(10,5))
        grid = [fig.add_subplot(111)]
    else:
        fig, grid = plt.subplots(nrows=size, ncols=1, figsize=(10,6),
                                 sharex=True)
    i=0
    for P in predictors:
        ax=grid[i]
        P.plotBinders(name, ax=ax, n=n, **kwargs)
        i+=1
    ax.xaxis.set_major_locator(MaxNLocator(20))
    plt.tight_layout()
    print fig
    if save == True:
        filename = os.path.join(path,'%s.png' %name)
        fig.savefig(filename, dpi=dpi)
        plt.close(fig)
    else:
        filename = ''
    return fig, filename

def splitPredictionsFile(infile):
    """Split predictions into one file per protein"""

    path = os.path.splitext(infile)[0]
    if not os.path.exists(path):
        os.mkdir(path)
    res = pd.read_msgpack(infile,iterator=True)
    print 'read %s' %infile
    for df in res:
        if type(df) is not pd.DataFrame:
           continue
        name = df.name.values[0]
        outfile = os.path.join(path, name+'.mpk')
        print outfile
        pd.to_msgpack(outfile,df)
    return

def getNearest(df):
    """Get nearest binder"""
    grps = df.groupby('name')
    new = []
    def closest(x,g):
        if len(g.pos)==1:
            return 1
        return min([abs(x.pos-i) for i in g.pos if i!=x.pos])
    for i,g in grps:
        positions = g.pos
        g['nearest'] = g.apply(lambda x: closest(x,g),axis=1)
        new.append(g)
    df = pd.concat(new)
    return df

def getBinders(preds,n=3):
    """Get binders for multiple predictors"""
    b={}
    for m in preds:
        pred = preds[m]
        binders = pred.getPromiscuousBinders(n=n)
        if len(binders)>0:
            binders = binders.sort('pos')
            b[m] = binders
    return b

def getStandardmhc1Name(name):
    """Taken from iedb mhc1 utils.py"""
    temp = name.strip().split('-')
    length = temp[-1]
    mhc = '-'.join(temp[0:-1])
    return mhc

class Predictor(object):
    """Base class to handle generic predictor methods, usually these will
       wrap methods from other modules and/or call command line predictors.
       Subclass for specific functionality"""
    def __init__(self, data=None):
        self.data = data
        self.name = ''
        self.scorekey = 'score'
        self.operator = '<'
        self.rankascending = 1
        #can specify per allele cutoffs here
        self.allelecutoffs = {}
        return

    def predict(self, sequence, peptide):
        """Does the actual scoring. Must override this.
           Should return a pandas DataFrame"""
        return

    def prepareData(self, result, name, allele):
        """Put raw prediction data into DataFrame and rank,
           override for custom processing"""
        df = pd.DataFrame(result, columns=['peptide','core','pos','score'])
        df['name'] = name
        df['allele'] = allele
        self.getRanking(df)
        return df

    def getRanking(self, df):
        """Add a ranking column according to scorekey"""
        s=self.scorekey
        df['rank'] = df[s].rank(method='min',ascending=self.rankascending)
        df.sort_index(by=['rank','name','allele'], ascending=True, inplace=True)
        return

    def evaluate(self, df, key, value, operator='<'):
        """Evaluate binders depending on less than or greater than a
           cutoff"""
        if operator == '<':
            return df[df[key] <= value]
        else:
            return df[df[key] >= value]

    def getBinders(self, method='cutoff', q=0.01, data=None):
        """Return top % ranked or using cutoff"""
        if not data is None:
            df = data
        else:
            df = self.data
        #print df
        key = self.scorekey
        op = self.operator
        if method == 'cutoff':
            #evaluate per allele in case we have tailored cutoffs
            res = []
            for a,g in df.groupby('allele'):
                if self.allelecutoffs.has_key(a):
                    cutoff = self.allelecutoffs[a]
                else:
                    cutoff = self.cutoff
                b = self.evaluate(g, key, cutoff, op)
                res.append(b)
            return pd.concat(res)
            #return self.evaluate(df, key,self.cutoff,op)
        elif method == 'rank':
            #get top ranked % per protein
            res=[]
            for i,g in df.groupby('name'):
                value = g['rank'].quantile(q=q)
                b = g[g['rank'] <= value]
                res.append(b)
            return pd.concat(res)

    def getPromiscuousBinders(self, n=3, method='cutoff', data=None, name=None):
        """Return only top binders present in at least n alleles"""

        if data is None:
            data = self.data
        else:
            self.data = data
        if name != None:
            data = self.data[self.data.name==name]
        df = self.getBinders(method, data=data)
        grps = df.groupby(['peptide','pos','name'])
        if self.operator == '<':
            func = min
        else:
            func = max
        s = grps.agg({'allele':pd.Series.count,self.scorekey:func})
        s = s[s.allele>=n]
        s = s.reset_index()
        #merge frequent binders with original data to retain fields
        p = list(data.groupby('allele'))[0][1]
        p = p.drop(['allele','rank',self.scorekey],1)

        if not s.empty:
            final = pd.merge(p,s,how='right',on=['peptide','pos','name'])
            l = self.getLength()
            #if l > 9:
            g = final.groupby('core')
            final = g.agg({self.scorekey:max,'name':first,'peptide':first,'pos':first,
                            'allele':first})
            final = final.reset_index().sort('pos')
            #name = df.head(1).name.max()
            #print merged.sort('pos')
            #print final
            return final
        #if not s.empty:
        #    return pd.merge(p,s,how='right',on=['peptide','pos','name'])
        else:
            return pd.DataFrame()

    def getUniqueCores(self, binders=False):
        """Get only unique cores"""
        if binders == True:
            df = self.getBinders()
        else:
            df = self.data
        grouped = df.groupby('core')
        cores = grouped.agg({self.scorekey:max})
        #cores = df.loc[grouped[self.scorekey].max().index]
        cores.sort(self.scorekey, inplace=True, ascending=self.rankascending)
        #print cores
        return cores

    def predictSequences(self, data, seqkey='peptide', length=11,
                        alleles=['HLA-DRB1*0101'], save=False):
        results=[]
        for i,row in data.iterrows():
            seq = row[seqkey]
            if len(seq)<=length: continue
            #print i,seq
            res=[]
            for a in alleles:
               df = self.predict(sequence=seq,length=length,
                                    allele=a,name=i)
               res.append(df)
            res = pd.concat(res)
            results.append(res)
            if save==True:
                pd.to_msgpack('predictions_%s.mpk' %self.name, res, append=True)
        self.data = pd.concat(results)
        return results

    def predictProteins(self, recs, length=11, names=None,
                         alleles=[], save=False, label='', path=''):
        """Get predictions for a set of proteins and/or over multiple alleles
          recs: a pandas DataFrame with cds data
          returns a dataframe of predictions over multiple proteins"""

        if type(alleles) is not types.ListType:
            alleles = [alleles]
        self.length = length
        if save == True:
            fname = 'epit_%s_%s_%s.mpk' %(label,self.name,length)
            fname = os.path.join(path,fname)
            print 'results will be saved to %s' %os.path.abspath(fname)
            meta = {'method':self.name, 'length':self.length}
            pd.to_msgpack(fname, meta)
        recs = Genome.getCDS(recs)
        if names != None:
            recs = recs[recs.locus_tag.isin(names)]
        proteins = list(recs.iterrows())
        results=[]
        for i,row in proteins:
            st=time.time()
            seq = row['translation']
            name = row['locus_tag']
            print name
            res = []
            for a in alleles:
                df = self.predict(sequence=seq,length=length,
                                    allele=a,name=name)
                if df is not None:
                    res.append(df)
            res = pd.concat(res)
            if save == True:
                pd.to_msgpack(fname, res, append=True)
            if save == False and len(df)>0:
                results.append(res)
        if save == False and len(results)>0:
            self.data = pd.concat(results)
        return

    def save(self, label, singlefile=True):
        """Save all current predictions dataframe with some metadata"""

        if singlefile == True:
            fname = 'epit_%s_%s_%s.mpk' %(label,self.name,self.length)
            print 'saving as %s' %fname
            meta = {'method':self.name, 'length':self.length}
            pd.to_msgpack(fname, meta)
            for i,g in self.data.groupby('name'):
                pd.to_msgpack(fname, g, append=True)
        else:
            #save one file per protein/name
            path = os.path.join(label,self.name)
            print 'saving to %s' %path
            if not os.path.exists(path):
                os.makedirs(path)
            for name,df in self.data.groupby('name'):
                outfile = os.path.join(path, name+'.mpk')
                pd.to_msgpack(outfile,df)
        return

    def getLength(self):
        """Get peptide length of current set of predictions"""
        if len(self.data)>0:
            return len(self.data.head(1).peptide.max())
        return

    def summary(self):
        '''print 'high binders: %s' %len(self.getBinders())
        print 'binders with unique cores: %s' %len(self.getUniqueCores(binders=True))
        allelegrps = self.data.groupby('allele')
        print '%s peptides in %s proteins and %s alleles' %(len(self.data),
                                            len(proteins),len(allelegrps))'''
        return

    def reshape(self, name=None):
        """Return pivoted data over alleles for summary use"""
        df = self.data
        if name != None:
            df = df[df.name==name]
        p = df.pivot(index='peptide', columns='allele', values=self.scorekey)
        p = p.reset_index()
        x = list(df.groupby('allele'))[0][1]
        p = p.merge(x[['pos','peptide']],on='peptide')
        p['mean'] = p.mean(1)
        p=p.sort('mean',ascending=False)
        return p

    def getNames(self):
        grp = self.data.groupby('name')
        return sorted(dict(list(grp)).keys())

    def benchmark(self):
        """Benchmark on known cores"""
        hits=0; total=0
        templates = Threading.templates
        for allele,row in templates.iterrows():
            name = row['pdbid']
            nativecore = row['core']
            seq = row['peptide']
            if len(seq)<=11: continue
            df = self.predict(seq,allele=allele,length=9)
            if len(df)==0: continue
            rank = df[df.peptide==nativecore]['rank'].values[0]
            #print df
            print allele,df.iloc[0].peptide,nativecore,rank,df.iloc[0][self.scorekey]
            if rank==1:
                hits+=1
            total+=1
        print '%s/%s correct' %(hits,total)
        return

    def plotBinders(self, name, cldist=7, n=2, tmregions=None,
                    legend=False, figsize=(9,3), ax=None):
        """Plot binders as bars per allele"""

        fig=None
        if ax==None:
            fig=plt.figure(figsize=figsize)
            ax=fig.add_subplot(111)
        cmap = plt.cm.get_cmap('jet')
        sckey = self.scorekey
        df = self.data[self.data.name==name]
        if self.cutoff < 0:
            highest = min(df[sckey])
        else:
            highest = max(df[sckey])
        if len(self.allelecutoffs)>0:
            lowest = min(self.allelecutoffs.values())
        else:
            lowest = self.cutoff
        lims = (lowest,highest)
        pb = self.getPromiscuousBinders(data=df,n=n)
        #pball = self.getPromiscuousBinders(data=df,n=1)
        grps = df.groupby('allele')
        cl = dbscan(pb,dist=cldist)
        j=0
        labels = []
        leg = []
        if len(pb)>0:
            for a, df in grps:
                c = cmap(float(j)/(len(grps)))
                b = self.getBinders(data=df)
                ind = np.arange(len(df))
                b = b[b.pos.isin(pb.pos)] #show only promiscuous
                b.sort('pos',inplace=True)
                y = b[sckey].values
                x = b['pos'].values
                bars = plotBars(x,y,ind,color=c,ax=ax,label='')
                labels.extend(zip(bars,b.peptide))
                if len(bars)>0:
                    leg.append((a,bars[0]))
                j+=1
        ax.set_title(self.name+' '+name)
        ax.set_ylim(lims)
        plt.setp(ax.get_xticklabels(), visible=True)

        #moving average plot of epitope density
        #df=df.sort('pos')
        #m = df[sckey].where(df[sckey].isin(pb.pos),other=0)
        #y = m.apply(lambda x: pball.sco
        #y = pd.stats.moments.rolling_mean(m, 10)
        #ax2 = ax.twinx()
        #ax2.plot(df.pos.values, y, '-',lw=2)

        plotRegions(cl,ax,alpha=0.2,ls='dashed')
        if tmregions != None:
             plotRegions(tmregions,ax,color='y',alpha=0.2,ls='dashed')
        if legend == True and len(leg)>0:
            patches,l = zip(*leg)
            ax.legend(l,patches,fontsize=9,mode="expand",ncol=6,framealpha=0.5)
        plt.tight_layout()
        return fig, labels

class NetMHCIIPanPredictor(Predictor):
    """netMHCIIpan predictor"""
    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'netmhciipan'
        self.colnames = ['pos','HLA','peptide','Identity','Pos','Core',
                         '1-log50k(aff)','Affinity','Rank']
        self.scorekey = '1-log50k(aff)'
        self.cutoff = 0.2
        self.operator = '>'
        self.rankascending = 0

    def readResult(self, res):
        """Read results from netMHCIIpan"""
        data=[]
        res = res.split('\n')[43:]
        ignore=['Protein','pos','']
        for r in res:
            if r.startswith('-'): continue
            row = re.split('\s*',r.strip())[:9]
            if len(row)!=9 or row[0] in ignore:
                continue
            data.append(dict(zip(self.colnames,row)))
        return data

    def prepareData(self, df, name):

        df = df.convert_objects(convert_numeric=True)
        df['name'] = name
        df.rename(columns={'Core': 'core','HLA':'allele'}, inplace=True)
        df=df.drop(['Pos','Identity','Rank'],1)
        self.getRanking(df)
        self.data = df
        return

    def runSequence(self, seq, length, allele):
        seqfile = createTempSeqfile(seq)
        cmd = 'netMHCIIpan -s -l %s -a %s %s' %(length, allele, seqfile)
        temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        rows = self.readResult(temp)
        df = pd.DataFrame(rows)
        return df

    def predict(self, sequence=None, peptides=None, length=11,
                    allele='HLA-DRB1*01:01', name='',
                    pseudosequence=None):
        """Call netMHCIIpan command line"""

        allele = allele.split('-')[1].replace('*','_')
        if peptides != None:
            res = pd.DataFrame()
            for p in peptides:
                temp = self.runSequence(p, len(p), allele)
                res = res.append(temp,ignore_index=True)
        else:
            res = self.runSequence(sequence, length, allele)
        if len(res)==0:
            return res
        self.prepareData(res, name)
        #print self.data[self.data.columns[:7]][:5]
        return self.data

class IEDBMHCIPredictor(Predictor):
    """Using IEDB tools method, requires iedb-mhc1 tools"""
    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'iedbmhc1'
        self.scorekey = 'percentile_rank'
        self.cutoff = 3.0
        self.operator = '<'
        self.rankascending = 1
        self.methods = ['ann','arb','comblib_sidney2008','consensus',
                        'IEDB_recommended','netmhcpan','smm']
        self.path = '/local/iedbmhc1/'

    def predict(self, sequence=None, peptides=None, length=11,
                   allele='HLA-A*01:01', method='IEDB_recommended', name=''):
        """Use iedb MHCII python module to get predictions.
           Requires that the iedb MHC tools are installed locally"""

        seqfile = createTempSeqfile(sequence)
        cmd = os.path.join(self.path,'src/predict_binding.py')
        cmd = cmd+' %s %s %s %s' %(method,allele,length,seqfile)
        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except Exception as e:
            return None
        self.prepareData(temp, name)
        return self.data

    def prepareData(self, rows, name):

        df = pd.read_csv(StringIO.StringIO(rows),sep="\t")
        extracols = ['ann_ic50','ann_rank','smm_ic50','smm_rank',
                    'comblib_sidney2008_score','comblib_sidney2008_rank',
                    'netmhcpan_rank','percentile_rank']
        df = df.drop(extracols,1)
        df.reset_index(inplace=True)
        df.rename(columns={'index':'pos',
                           'percentile_rank':'method',
                           'method':'percentile_rank'},
                           inplace=True)
        df['core'] = df.peptide
        df['name'] = name
        df.sort_index(by=['percentile_rank'], ascending=True, inplace=True)
        self.getRanking(df)
        self.data = df
        #print self.data[:10]
        return

    def getMHCIList(self):
        """Get available alleles from model_list file and
            convert to standard names"""
        afile = os.path.join(self.path,'data/MHCI_mhcibinding20130222/consensus/model_list.txt')
        df = pd.read_csv(afile,sep='\t',names=['name','x'])
        alleles = list(df['name'])
        alleles = sorted(list(set([getStandardmhc1Name(i) for i in alleles])))
        return alleles

class IEDBMHCIIPredictor(Predictor):
    """Using IEDB mhcii method, requires iedb-mhc2 tools"""
    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'iedbmhc2'
        self.scorekey = 'consensus_percentile'
        self.cutoff = 3
        self.operator = '<'
        self.rankascending = 1
        self.methods = ['arbpython','comblib','consensus3','IEDB_recommended',
                    'NetMHCIIpan','nn_align','smm_align','tepitope']
        self.path = '/local/iedbmhc2/'

    def prepareData(self, rows, name):
        df = pd.read_csv(StringIO.StringIO(rows),delimiter=r"\t")
        extracols = ['Start','End','comblib_percentile','smm_percentile','nn_percentile',
                'Sturniolo core',' Sturniolo score',' Sturniolo percentile']
        df = df.drop(extracols,1)
        df.reset_index(inplace=True)
        df.rename(columns={'index':'pos','Sequence': 'peptide','Allele':'allele'},
                           inplace=True)
        df['core'] = df.nn_core
        df['name'] = name
        self.getRanking(df)
        self.data = df
        return

    def predict(self, sequence=None, peptides=None, length=15,
                   allele='HLA-DRB1*01:01', method='consensus3', name=''):
        """Use iedb MHCII python module to get predictions.
           Requires that the iedb MHC tools are installed locally"""

        seqfile = createTempSeqfile(sequence)
        cmd = os.path.join(self.path,'mhc_II_binding.py')
        cmd = cmd+' %s %s %s' %(method,allele,seqfile)
        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except:
            print 'allele %s not available?' %allele
            return None
        self.prepareData(temp, name)
        #print self.data
        return self.data

class TEpitopePredictor(Predictor):
    """Predictor using Tepitope QM method"""
    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'tepitope'
        self.pssms = Tepitope.getPSSMs()
        self.cutoff = 2
        self.operator = '>'
        self.rankascending = 0

    def predict(self, sequence=None, peptides=None, length=11,
                    allele='HLA-DRB1*0101', name='',
                    pseudosequence=None):

        self.sequence = sequence
        if not allele in self.pssms:
            #print 'computing virtual matrix for %s' %allele
            try:
                m = Tepitope.createVirtualPSSM(allele)
            except:
                return pd.DataFrame()
        else:
            m = self.pssms[allele]
        m = m.transpose().to_dict()
        result = Tepitope.getScores(m, sequence, peptides, length)
        df = self.prepareData(result, name, allele)
        self.data = df
        #print df[:12]
        return df

class ThreadingPredictor(Predictor):
    """Predictor using basic threading, templates should have peptides
       with length 9,11,13 or 15"""
    def __init__(self, data=None, template=None):
        Predictor.__init__(self, data=data)
        self.name = 'threading'
        self.cutoff = -5
        self.templatecontacts = {}
        self.rankascending = 1

    def predictProteins(self, recs, **kwargs):
        #pre-calculate the template interactions for each allele
        alleles = kwargs['alleles']
        self.templatecontacts = Threading.preCalculateContacts(alleles)
        Predictor.predictProteins(self, recs=recs, **kwargs)
        return

    def predict(self, sequence=None, peptides=None, length=9,
                 allele='HLA-DRB1*0101', name=''):
        """Predict threading scores using available strutural templates"""

        self.sequence = sequence
        #print allele
        tc = self.templatecontacts
        if allele not in tc:
            template = Threading.getTemplateStructure(allele)
            if template == None:
                template = Threading.modelAlleleStructure(allele)

            inter = Threading.getPocketContacts(template)
            tc[allele] = inter
        else:
            inter = tc[allele]
        result = Threading.getScores(sequence, peptides,
                                     inter=inter, length=length)
        df = self.prepareData(result, name, allele)
        self.data = df
        #print df
        return self.data


