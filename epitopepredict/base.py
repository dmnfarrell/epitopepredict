#!/usr/bin/env python

"""
    MHC prediction base module for core classes
    Created November 2013
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, shutil, string
import csv, glob, pickle
import time, io
import operator as op
import re, types
import math
import subprocess
from subprocess import CalledProcessError
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from . import utilities, peptutils, sequtils, tepitope

home = os.path.expanduser("~")
path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(path, 'mhcdata')
predictors = ['tepitope','netmhciipan','iedbmhc1','iedbmhc2','bcell']
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

#these paths should be set by user before calling predictors
iedbmhc1path = '/local/iedbmhc1/'
iedbmhc2path = '/local/iedbmhc2/'
iedbbcellpath = '/local/iedbbcell/'

#six class I super-type alleles
mhc1supertypes = ['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*03:01', 'HLA-A*24:02',
                  'HLA-B*07:02','HLA-B*44:03']
#eight class II super-types
mhc2supertypes = ['HLA-DRB1*0101','HLA-DRB1*0301','HLA-DRB1*0401','HLA-DRB1*0701'
                  'HLA-DRB1*0801','HLA-DRB1*1101','HLA-DRB1*1301','HLA-DRB1*1501']

def first(x):
    return x.iloc[0]

def getIEDBRequest(seq, alleles='HLA-DRB1*01:01', method='consensus3'):
    import requests
    url = 'http://tools.iedb.org/tools_api/mhcii/'
    values = {'method' : method,
              'sequence_text' : seq,
              'allele' : alleles }
    r=requests.post(url,data=values)
    df=pd.read_csv(io.StringIO(r.content),sep='\t')
    #df=df.drop(['nn_align_core','nn_align_ic50','nn_align_rank'])
    return df

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
    """Get a predictor"""

    if name == 'netmhciipan':
        return NetMHCIIPanPredictor(**kwargs)
    elif name == 'iedbmhc1':
        return IEDBMHCIPredictor(**kwargs)
    elif name == 'iedbmhc2':
        return IEDBMHCIIPredictor(**kwargs)
    elif name == 'bcell':
        return BCellPredictor(**kwargs)
    elif name == 'tepitope':
        return TEpitopePredictor(**kwargs)
    elif name == 'mhcflurry':
        return MHCFlurryPredictor(**kwargs)
    else:
        print ('no such predictor %s' %name)
        return

def getLength(data):
    """Get peptide length of a dataframe of predictions"""

    if len(data)>0:
        return len(data.head(1).peptide.max())
    return

def getCoordsfromSequence(df, genome):
    """Get peptide coords from parent protein sequences"""

    def func(x):
        seq = x.translation
        st = seq.find(x.peptide)
        end = st+len(x.peptide)
        return pd.Series({'start':st,'end':end})
    temp = df.merge(genome[['locus_tag','translation']],
                    left_on='name',right_on='locus_tag')
    temp = temp.apply( lambda r: func(r),1)
    df = df.drop(['start','end'],1)
    return df.join(temp)

def getCoords(df):
    """Get start end coords from position and length of peptides"""

    if 'start' in df.columns:
        return df
    df['start'] = df.pos.astype(int)
    df['end'] = ( df.pos + df.peptide.str.len() ).astype(int)
    return df

def createTempSeqfile(sequences, seqfile='tempseq.fa'):

    if isinstance(sequences, str):
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

def getSequencefromPredictionData(data):
    """Guess original sequence from predictions dataframe"""

    return seq

'''def getOverlappingBinders(B):
    for df in B:
        df.set_index('peptide',inplace=True)
    print (B)
    x = pd.concat(B,join='inner')
    #x =  x[['pos']].sort('pos')
    return x'''

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

'''def getBinders(preds,n=3):
    """Get binders for multiple predictors"""

    b={}
    for m in preds:
        pred = preds[m]
        binders = pred.getPromiscuousBinders(n=n)
        if len(binders)>0:
            binders = binders.sort('pos')
            b[m] = binders
    return b'''

def summarize(data):
    """Summarise prediction data"""

    #print 'binders with unique cores: %s' %len(self.getUniqueCores(binders=True))
    allelegrps = data.groupby('allele')
    proteins = data.groupby('name')
    print ('summary: %s peptides in %s proteins and %s alleles' %(len(data),
                                        len(proteins),len(allelegrps)))
    return

def getBindersfromPath(method, path, n=3, cutoff_method='default',
                       perc=0.98, promiscuous=True):
    """
    Get all binders from a set of binding results stored in a directory.

    Args:
        path: The file path with all the binding prediction results
        cutoff_method: Prediction method used to create the data
        n: minimum number of alleles if using promiscuous binders
        perc: percentile for cutoff(s)
        promiscuous: whether to return only promiscuous binders

    Returns:
        A dataframe with all binders matching the required critieria
    """

    binders = []
    P = getPredictor(method)
    files = glob.glob(os.path.join(path, '*.csv'))

    #get allele specific cutoffs
    if cutoff_method == 'default':
        P.allelecutoffs = getCutoffs(P, path, perc, overwrite=True)

    key = P.scorekey
    for f in files:
        df = pd.read_csv(f, index_col=0)
        if not key in df.columns:
            continue
        #print (df[:3])
        if promiscuous== True:
            b = P.getPromiscuousBinders(data=df, n=n,
                                        cutoff_method=cutoff_method)
        else:
            b = P.getBinders(data=df)
        #print b[:5]
        binders.append(b)
    result = pd.concat(binders)
    return result

def getCutoffs(predictor, path=None, perc=0.98, overwrite=True):
    """
    Estimate global allele-based cutoffs for predictions.
    Args:
        predictor: A Predictor object with data loaded optionally
        path: An optional file path with the binding prediction results
        method: Prediction method
        perc: percentile level of score to select cutoffs
    Returns:
        A dictionary with cutoff values per allele
    """

    if path != None:
        binsfile = os.path.join(path, 'quantiles.csv')
    else:
        binsfile = ''
    if not os.path.exists(binsfile) or overwrite==True:
        bins = getScoreDistributions(predictor, path)
        if bins is None:
            return {}
        if binsfile != '':
            bins.to_csv(binsfile, float_format='%.3f')
    else:
        bins = pd.read_csv(binsfile, index_col=0)
    cutoffs = dict(bins.ix[perc])
    return cutoffs

def getScoreDistributions(predictor, path=None):
    """Get global score distributions and save quantile values for each allele
       Assumes all the files in path represent related proteins"""

    if path != None:
        predictor.load(path=path, file_limit=500)
    df = predictor.data
    if df is None or len(df)==0:
        print ('no prediction data loaded')
        return
    key = predictor.scorekey
    x = df.pivot_table(index='peptide', columns='allele', values=key)
    percs = np.arange(0.01,1,0.01)
    bins = x.quantile(percs)
    #reverse if best values are lower
    if predictor.operator == '<':
        bins.index = pd.Series(bins.index).apply(lambda x: 1-x)
    return bins

def getStandardMHCI(name):
    """Taken from iedb mhc1 utils.py"""

    temp = name.strip().split('-')
    length = temp[-1]
    mhc = '-'.join(temp[0:-1])
    return mhc

def getDRBList(a):
    """Get DRB list in standard format"""

    s = pd.Series(a)
    s = s[s.str.contains('DRB')]
    s = s.apply(lambda x:'HLA-'+x.replace('_','*'))
    return list(s)

def getDQPList(a):
    """Get DRB list in standard format"""
    s = pd.Series(a)
    s = s[s.str.contains('DQ')]
    s = s.apply(lambda x:x.replace('_','*'))
    return list(s)

def getStandardMHCII(x):
    return 'HLA'+x.replace('_','*')

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

    def __repr__(self):

        if (self.data is None) or len(self.data) == 0:
            return '%s predictor' %self.name
        else:
            n = len(self.data.name.unique())
            return '%s predictor with results in %s proteins' %(self.name, n)

    def predict(self, sequence, peptide):
        """Does the actual scoring. Must override this.
           Should return a pandas DataFrame"""
        return

    def prepareData(self, result, name, allele):
        """Put raw prediction data into DataFrame and rank,
           override for custom processing. Can be overriden for
           custom data."""

        df = pd.DataFrame(result, columns=['peptide','core','pos','score'])
        df['name'] = name
        df['allele'] = allele
        self.getRanking(df)
        return df

    def getRanking(self, df):
        """Add a ranking column according to scorekey"""

        s=self.scorekey
        df['rank'] = df[s].rank(method='min',ascending=self.rankascending)
        df.sort_values(by=['rank','name','allele'], ascending=True, inplace=True)
        return

    def evaluate(self, df, key, value, operator='<'):
        """
        Evaluate binders less than or greater than a cutoff.
        This method is called by all predictors to get binders
        """

        if operator == '<':
            return df[df[key] <= value]
        else:
            return df[df[key] >= value]

    def getBinders(self, by='score', cutoff_method='default', perc=0.98,
                   data=None, name=None):
        """
        Get the top scoring binders using percentile ranking or single cutoff.
        Args:
            data: binding predictions for one or more proteins
            cutoff_method: method to use for binder thresholds
            perc: percentile threshold for ranking in each allele
        Returns:
            pandas DataFrame of all binders
        """

        if data is None:
            if self.data is None:
                print ('no prediction data available')
                return
            else:
                data = self.data
        if name != None:
            if name not in list(data.name):
                print ('no such protein in data')
                return
            data = data[data.name==name]

        key = self.scorekey
        op = self.operator
        if op == '<':
            q = 1-perc
        else:
            q = perc

        cutoffs = {}
        if cutoff_method == 'global':
            #we derive cutoffs using all loaded data
            for a,g in self.data.groupby('allele'):
                cutoffs[a] = g[key].quantile(q=q)
        else:
            cutoffs = self.allelecutoffs
        if cutoff_method == 'simple':
            #we just use a single cutoff value for all
            cutoff = self.cutoff
            b = self.evaluate(data, key, cutoff, op)
            return b
        else:
            res=[]
            for a,g in data.groupby('allele'):
                if a in cutoffs:
                    cutoff = cutoffs[a]
                else:
                    cutoff = g[key].quantile(q=q)
                #print (a,perc,value)
                b = self.evaluate(g, key, cutoff, op)
                if b is not None:
                    res.append(b)
            if len(res) > 0:
                return pd.concat(res)
        return

    def getPromiscuousBinders(self, n=2, cutoff_method='default', perc=0.98,
                              data=None, name=None):
        """
        Get top scoring binders present in at least n alleles.
        Args:
            n: number of alleles
            cutoff method: method to use for cutoffs - default or simple
            data: a dataframe of prediction data, optional
            name: name of the proteins to use, required if data contains multiple proteins
            perc: percentile cutoff, applied per allele
        Returns:
            pandas DataFrame with binders
        """

        if data is None:
            data = self.data
        #get binders using the provided or current prediction data
        b = self.getBinders(cutoff_method, data=data, name=name, perc=perc)
        if b is None or len(b) == 0:
            return pd.DataFrame()
        grps = b.groupby(['peptide','pos','name'])
        if self.operator == '<':
            func = min
            skname = 'min'
        else:
            func = max
            skname = 'max'
        s = grps.agg({'allele':pd.Series.count, self.scorekey:[func,np.mean]})
        s.columns = s.columns.get_level_values(1)
        s.rename(columns={skname: self.scorekey, 'count': 'alleles'}, inplace=True)

        s = s[s.alleles>=n]
        s = s.reset_index()
        #merge frequent binders with original data to retain fields
        p = list(data.groupby('allele'))[0][1]
        p = p.drop(['allele','rank',self.scorekey],1)

        if not s.empty:
            final = pd.merge(p,s,how='right',on=['peptide','pos','name'])
            l = getLength(b)
            grpkey = 'core'
            if not grpkey in final.columns:
                final[grpkey] = final.peptide
            g = final.groupby(grpkey)
            final = g.agg({self.scorekey:max, 'name':first, 'peptide': first,
                        'pos':first, 'alleles':first, 'mean':first})
            final = final.reset_index().sort_values('pos')
            #print final
            return final
        else:
            return pd.DataFrame()

    def consensusRankedBinders(self, name=None, how='median', threshold=None):
        """
        Get the top percentile rank of each binder over all alleles.
        Args:
            name: specify protein name, otherwise all current data used
            how: method to use for rank selection, 'median', 'best' or
                 'mean'
            threshold: apply a threshold
        """

        df = self.data
        if name != None:
            df=df[df.name==name]
        funcs = { 'median':np.median, 'mean':np.mean, 'best':min }
        func = funcs[how]
        b = df.groupby(['peptide']).agg({'rank': func,'pos':first, 'name':first})
        b = b.reset_index().sort_values('rank')
        if threshold != None:
            b = b[b['rank'] < threshold]
        return b

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

    def predictPeptides(self, peptides, length=11,
                        alleles=[], save=False):
        """Predict a set of individual peptides"""

        results=[]
        for seq in peptides:
            print (seq)
            if len(seq)<length: continue
            res=[]
            for a in alleles:
               df = self.predict(sequence=seq, length=length,
                                    allele=a, name=seq)
               res.append(df)
            res = pd.concat(res)
            results.append(res)
        self.data = pd.concat(results)
        return results

    def predictProteins(self, recs, key='locus_tag', seqkey='translation',
                        length=11, overlap=1, names=None, alleles=[],
                        path=None, overwrite=True, cpu=1):
        """Get predictions for a set of proteins and/or over multiple alleles
          Args:
            recs: protein sequences in a pandas DataFrame
            length: length of peptides to predict
            names: names of proteins to use from sequences
            alleles: allele list
            path: if results are to be saved to disk provide a path, otherwise results
            for all proteins are stored in the data attribute of the predictor
            overwrite: over write existing protein files in path if present
            cpu: number of processors to use
          Returns:
            a dataframe of predictions over multiple proteins"""

        if type(alleles) is str:
            alleles = [alleles]
        elif type(alleles) is pd.Series:
            alleles = alleles.tolist()
        if len(alleles) == 0:
            print ('no alleles provided')
            return
        self.length = length
        #recs = sequtils.getCDS(recs)

        if names != None:
            recs = recs[recs[seqkey].isin(names)]
        proteins = list(recs.iterrows())
        results = []
        if path is not None and path != '':
            if not os.path.exists(path):
                os.mkdir(path)

        if cpu > 1:
            #use more than one cpu
            from multiprocessing import Process,Queue
            grps = np.array_split(proteins, cpu)
            procs = []
            queue = Queue()
            for df in grps:
                p = Process(target=self._predict_multiple,
                            kwargs={'proteins':df, 'alleles':alleles,
                                    'path': path, 'overwrite':overwrite,
                                    'length':length, key:key, seqkey:seqkey,
                                    'queue':queue})
                procs.append(p)
                p.start()
            #print ( [queue.get() for p in procs])
            for p in procs:
                p.join()
            #print (len(results))
        else:
            results = self._predict_multiple(proteins, path, overwrite, alleles, length,
                                             overlap=overlap, key=key, seqkey=seqkey)

        print ('predictions done for %s proteins in %s alleles' %(len(proteins),len(alleles)))
        if path is None:
            #if no path we keep assign results to the data object
            #assumes we have enough memory..
            self.data = results
        else:
            print ('results saved to %s' %os.path.abspath(path))
        return

    def _predict_multiple(self, proteins, path, overwrite, alleles, length, overlap=1,
                          key='locus_tag', seqkey='sequence', queue=None):
        """Predictions for multiple proteins in a dataframe
            Args: as for predictProteins
            Returns: a dataframe of the results if no path is given
        """

        results = []
        for i,row in proteins:
            seq = row[seqkey]
            name = row[key]
            if path is not None:
                fname = os.path.join(path, name+'.csv')
                if os.path.exists(fname) and overwrite == False:
                    continue
            res = []
            for a in alleles:
                df = self.predict(sequence=seq,length=length,overlap=overlap,
                                    allele=a,name=name)
                if df is not None:
                    res.append(df)
                #print (name, a, df.pos.max())
            res = pd.concat(res)
            #print (res.pos.max())
            if path is not None and len(res)>0:
                res.to_csv(fname)
            else:
                results.append(res)
        if len(results)>0:
            results = pd.concat(results)
        #if queue != None:
            #print (path,len(results))
            #queue.put(results)
        return results

    def load(self, filename=None, path=None, names=None,
               compression='infer', file_limit=None):
        """
        Load results for one or more proteins
        Args:
            filename: name of a csv file with predictions
            path: directory with one or more csv files
            file_limit: limit to load only the this number of proteins
        """

        if filename != None:
            self.data = pd.read_csv(filename, index_col=0)
        elif path != None:
            if not os.path.exists(path):
                print('no such path %s' %path)
                return
            files = glob.glob(os.path.join(path, '*.csv'))
            if names is not None:
                #if type(names) is pd.Series:
                #    names = list(names)
                names = [n+'.csv' for n in names]
                files = [n for n in files if os.path.basename(n) in names]
                #print (len(files))
            if file_limit != None:
                files = files[:file_limit]
            res = []
            for f in files:
                df = pd.read_csv(f, index_col=0, compression=compression)
                if len(df) == 0:
                    continue
                if not self.scorekey in df.columns:
                    continue
                res.append(df)
            self.data = pd.concat(res)
        return

    def save(self, prefix='_', filename=None, compression=None):
        """
        Save all current predictions dataframe with some metadata
        Args:
            prefix: if writing to a path, the prefix name
            filename: if saving all to a single file
            compression: a string representing the compression to use,
            allowed values are 'gzip', 'bz2', 'xz'.
        """

        exts = {'gzip':'.gz','bz2':'.bz2','xz':'.xz'}
        if filename != None:
            cext = exts[compression]
            if compression != None and not filename.endswith(cext):
                filename += cext
            self.data.to_csv(filename, compression=compression)
        else:
            #save one file per protein/name
            ext = '.csv'
            path = os.path.join(prefix, self.name)
            print ('saving to %s' %path)
            if not os.path.exists(path):
                os.makedirs(path)
            for name,df in self.data.groupby('name'):
                outfile = os.path.join(path, name+ext)
                df.to_csv(outfile)
        return

    def save_msgpack(self, filename=None):
        """Save as msgpack format - experimental"""

        if filename == None:
            filename = 'epit_%s_%s_%s.msg' %(label,self.name,self.length)
        print ('saving as %s' %filename)
        meta = {'method':self.name, 'length':self.length}
        pd.to_msgpack(filename, meta)
        for i,g in self.data.groupby('name'):
            pd.to_msgpack(filename, g, append=True)
        return

    def summarize(self):
        """Summarise currently loaded data"""

        return summarize(self.data)

    def alleleSummary(self, perc=0.98):
        """Allele based summary"""

        b = self.getBinders(perc=perc)
        summary = b.groupby('allele').agg({'peptide':np.size,self.scorekey:np.mean})
        return summary

    def reshape(self, name=None):
        """Return pivoted data over alleles for summary use"""

        df = self.data
        if name != None:
            df = df[df.name==name]
        p = df.pivot_table(index='peptide', columns='allele', values=self.scorekey)
        p = p.reset_index()
        x = list(df.groupby('allele'))[0][1]
        p = p.merge(x[['pos','peptide']],on='peptide')
        p['mean'] = p.mean(1)
        p=p.sort('mean',ascending=self.rankascending)
        return p

    def getNames(self):
        grp = self.data.groupby('name')
        return sorted(dict(list(grp)).keys())

    def plot(self, name, **kwargs):
        """
        Use module level plotting.mpl_plot_tracks method for predictor plot
        Args:
            name:
            n: min no. of alleles to be visible
            perc: percentile cutoff for score
            cutoff_method: method to use for cutoffs
        """

        from . import plotting
        if name == None:
            return
        plot = plotting.mpl_plot_tracks([self], name=name, **kwargs)
        return plot

    def getAlleles(self):
        """Get available alleles - override"""
        return []

class NetMHCIIPanPredictor(Predictor):
    """netMHCIIpan predictor"""

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'netmhciipan'
        self.colnames = ['pos','HLA','peptide','Identity','Pos','Core',
                         '1-log50k(aff)','Affinity','Rank']
        self.scorekey = '1-log50k(aff)'
        self.cutoff = .426
        self.operator = '>'
        self.rankascending = 0

    def readResult(self, res):
        """Read raw results from netMHCIIpan output"""

        data=[]
        res = res.split('\n')[19:]
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
        #df = df.apply(pd.to_numeric)#, errors='ignore')
        df['name'] = name
        df.rename(columns={'Core': 'core','HLA':'allele'}, inplace=True)
        df = df.drop(['Pos','Identity','Rank'],1)
        df = df.dropna()
        df['allele'] = df.allele.apply( lambda x: self.convert_allele_name(x) )
        self.getRanking(df)
        self.data = df
        return

    def runSequence(self, seq, length, allele, overlap=1):
        """Run netmhciipan for a single sequence"""

        seqfile = createTempSeqfile(seq)
        cmd = 'netMHCIIpan -s -length %s -a %s -f %s' %(length, allele, seqfile)
        #print cmd
        temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        rows = self.readResult(temp)
        df = pd.DataFrame(rows)
        return df

    def predict(self, sequence=None, peptides=None, length=11, overlap=1,
                    allele='HLA-DRB1*0101', name='',
                    pseudosequence=None):
        """Call netMHCIIpan command line"""

        #assume allele names are in standard format HLA-DRB1*0101
        try:
            allele = allele.split('-')[1].replace('*','_')
        except:
            print('invalid allele')
            return
        allele = allele.replace(':','')
        if peptides != None:
            res = pd.DataFrame()
            for p in peptides:
                temp = self.runSequence(p, len(p), allele)
                res = res.append(temp,ignore_index=True)
        else:
            res = self.runSequence(sequence, length, allele, overlap)
        if len(res)==0:
            return res
        self.prepareData(res, name)
        #print self.data[self.data.columns[:7]][:5]
        return self.data

    def getAlleles(self):
        """Get available alleles"""

        cmd = 'netMHCIIpan -list'
        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except:
            print('netmhciipan not installed?')
            return []
        alleles=temp.split('\n')[34:]
        #print sorted(list(set([getStandardmhc1Name(i) for i in alleles])))
        return alleles

    def convert_allele_name(self, r):
        """Convert allele names to internally used form"""

        if not r.startswith('HLA'):
            return 'HLA-'+r.replace(':','')
        else:
            return r.replace(':','')

class IEDBMHCIPredictor(Predictor):
    """Using IEDB tools method, requires iedb-mhc1 tools"""
    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'iedbmhc1'
        self.scorekey = 'score'
        self.methods = {'ANN':'ann_ic50','IEDB_recommended':'smm_ic50',
                         'Consensus (ANN,SMM)':'ann_ic50','NetMHCpan':'netmhcpan_ic50'}
        self.cutoff = .426
        self.operator = '>'
        self.rankascending = 0
        self.iedbmethod = 'IEDB_recommended'
        #self.path = iedbmhc1path

    def predict(self, sequence=None, peptides=None, length=11, overlap=1,
                   allele='HLA-A*01:01', name=''):
        """Use iedb MHCII python module to get predictions.
           Requires that the iedb MHC tools are installed locally"""

        seqfile = createTempSeqfile(sequence)
        path = iedbmhc1path
        if not os.path.exists(path):
            print ('iedb mhcI tools not found')
            return
        cmd = os.path.join(path,'src/predict_binding.py')
        cmd = cmd+' %s %s %s %s' %(self.iedbmethod,allele,length,seqfile)
        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash',
                stderr=subprocess.STDOUT)
        except CalledProcessError as e:
            print (e)
            return
        self.prepareData(temp, name)
        return self.data

    def prepareData(self, rows, name):
        """Prepare data from results"""

        df = pd.read_csv(io.BytesIO(rows),sep="\t")
        if len(df)==0:
            return
        df = df.replace('-',np.nan)
        df = df.dropna(axis=1,how='all')
        df.reset_index(inplace=True)
        df.rename(columns={'index':'pos',
                           'percentile_rank':'method',
                           'method':'percentile_rank'},
                           inplace=True)
        df['core'] = df.peptide
        df['name'] = name
        key = self.getScoreKey(df)
        df['ic50'] = df[key]
        df['score'] = df.ic50.apply( lambda x: 1-math.log(x, 50000))
        self.data = df
        self.getRanking(df)
        self.data = df
        return

    def getScoreKey(self, data):
        """Get iedbmhc1 score key from data"""

        m = data['method'].head(1).squeeze()
        key = self.methods[m]
        return key

    def getMHCIList(self):
        """Get available alleles from model_list file and
            convert to standard names"""

        afile = os.path.join(iedbmhc1path, 'data/MHCI_mhcibinding20130222/consensus/model_list.txt')
        df = pd.read_csv(afile,sep='\t',names=['name','x'])
        alleles = list(df['name'])
        alleles = sorted(list(set([getStandardMHCI(i) for i in alleles])))
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
        #self.path = '/local/iedbmhc2/'

    def prepareData(self, rows, name):
        df = pd.read_csv(io.StringIO(rows),delimiter=r"\t")
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
        path = iedbmhc2path
        if not os.path.exists(path):
            print ('iedb mhcII tools not found')
            return
        cmd = os.path.join(path,'mhc_II_binding.py')
        cmd = cmd+' %s %s %s' %(method,allele,seqfile)
        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except:
            print ('allele %s not available?' %allele)
            return
        self.prepareData(temp, name)
        #print self.data
        return self.data

class TEpitopePredictor(Predictor):
    """Predictor using tepitope QM method"""

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'tepitope'
        self.pssms = tepitope.getPSSMs()
        self.cutoff = 2
        self.operator = '>'
        self.rankascending = 0

    def predict(self, sequence=None, peptides=None, length=9, overlap=1,
                    allele='HLA-DRB1*0101', name='',
                    pseudosequence=None):

        self.sequence = sequence
        allele = allele.replace(':','')
        if not allele in self.pssms:
            #print 'computing virtual matrix for %s' %allele
            m = tepitope.createVirtualPSSM(allele)
            if m is None:
                print ('no such allele')
                return pd.DataFrame()
        else:
            m = self.pssms[allele]
        m = m.transpose().to_dict()
        result = tepitope.getScores(m, sequence, peptides, length, overlap=overlap)
        df = self.prepareData(result, name, allele)
        self.data = df
        #print df[:12]
        return df

    def getAlleles(self):
        return tepitope.getAlleles()

class BCellPredictor(Predictor):
    """Using IEDB tools methods, requires iedb bcell tools.
       see http://tools.immuneepitope.org/bcell """

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'iedbmhc1'
        self.scorekey = 'Score'
        self.methods = ['Chou-Fasman', 'Emini', 'Karplus-Schulz',
                        'Kolaskar-Tongaonkar', 'Parker', 'Bepipred']
        self.cutoff = 0.9
        self.operator = '>'
        self.rankascending = 0
        self.iedbmethod = 'Bepipred'
        self.path = iedbbcellpath

    def predict(self, sequence=None, peptides=None, window=None, name=''):
        """Uses code from iedb predict_binding.py """

        value = self.iedbmethod
        currpath=os.getcwd()
        os.chdir(self.path)
        sys.path.append(self.path)
        from src.BCell import BCell
        bcell = BCell()
        filepath = os.path.join(self.path,'bcell_scales.pickle')
        picklefile = open(filepath, 'rb')
        scale_dict = pickle.load(picklefile)
        bcell.scale_dict = scale_dict[value]
        if window==None:
            window = bcell.window
        center = "%d" %round(int(window)/2.0)
        if value == 'Emini':
            results = bcell.emini_method(value, sequence, window, center)
        elif value == 'Karplus-Schulz':
            results = bcell.karplusshulz_method(value, sequence, window, center)
        elif value == 'Kolaskar-Tongaonkar':
            results = bcell.kolaskartongaonkar_method(value, sequence, window, center)
        elif value == 'Bepipred':
            results = bcell.bepipred_method(value, sequence, window, center)
        else:
            results = bcell.classical_method(value, sequence, window, center)

        threshold = round(results[1][0], 3)
        temp=results[0]
        self.prepareData(temp, name)
        os.chdir(currpath)
        return self.data

    def prepareData(self, temp, name):

        df = pd.read_csv(temp,sep=",")
        if len(df)==0:
            return
        #df = df.replace('-',np.nan)
        df = df.dropna(axis=1,how='all')
        #df.reset_index(inplace=True)
        df['name'] = name
        self.data = df
        #print (df)
        return

    def predictProteins(self, recs, names=None, save=False,
                        label='', path='', **kwargs):
        """Get predictions for a set of proteins - no alleles so we override
        the base method for this too. """

        recs = sequtils.getCDS(recs)
        if names != None:
            recs = recs[recs.locus_tag.isin(names)]
        proteins = list(recs.iterrows())
        for i,row in proteins:
            seq = row['translation']
            name = row['locus_tag']
            #print (name)
            res = self.predict(sequence=seq,name=name)
            if save == True:
                #fname = os.path.join(path, name+'.mpk')
                #pd.to_msgpack(fname, res)
                fname = os.path.join(path, name+'.csv')
                res.to_csv(fname)

        return

class MHCFlurryPredictor(Predictor):
    """
    Predictor using MHCFlurry for MHC-I predictions. Requires you to
    install the python package and dependencies.
    see https://github.com/hammerlab/mhcflurry
    """

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'mhcflurry'
        self.cutoff = .426
        self.operator = '>'
        self.scorekey = 'score'
        self.rankascending = 0
        try:
            from mhcflurry import predict
        except:
            print ('mhcflurry not installed!')
        return

    def predict(self, sequence=None, peptides=None, length=11, overlap=1,
                      allele='HLA-A0101', name=''):
        """Uses mhcflurry python classes for prediction"""

        self.sequence = sequence
        from mhcflurry import predict
        if peptides == None:
            peptides, s = peptutils.createFragments(seq=sequence, length=length)
        scores=[]
        pos=0
        result = predict(alleles=[allele], peptides=peptides)
        df = self.prepareData(result, name, allele)
        self.data = df
        return df

    def prepareData(self, df, name, allele):
        """Post process dataframe to alter some column names"""

        df = df.rename(columns={'Allele':'allele','Peptide':'peptide'})
        df['name'] = name
        df['pos'] = df.index
        df['score'] = df['Prediction'].apply( lambda x: 1-math.log(x, 50000) )
        df['allele'] = df.allele.apply( lambda x: self.convert_allele_name(x) )
        self.getRanking(df)
        return df

    def convert_allele_name(self, r):
        return r[:5]+'*'+r[5:7]+':'+r[7:]

    def getAlleles(self):
        import mhcflurry
        return mhcflurry.Class1BindingPredictor.supported_alleles()
