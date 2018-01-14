#!/usr/bin/env python

"""
    MHC prediction base module for core classes
    Created November 2013
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, shutil, string
import csv, glob, pickle, tempfile
import time, io
import operator as op
import re, types
try:
    unicode = unicode
except NameError:
    #unicode is undefined, must be Python 3
    str = str
    unicode = str
    basestring = (str,bytes)
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
predictors = ['tepitope','netmhciipan','iedbmhc1','iedbmhc2','mhcflurry','iedbbcell']
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
iedbmhc1path = ''
iedbmhc2path = ''
iedbbcellpath = ''

mhc1_presets = ['mhc1_supertypes','us_caucasion_mhc1','us_african_mhc1','broad_coverage_mhc1']
mhc2_presets = ['mhc2_supertypes','human_common_mhc2','bovine_like_mhc2']

#sequence for testing
testsequence = ('MRRVILPTAPPEYMEAIYPVRSNSTIARGGNSNTGFLTPESVNGDTPSNPLRPIADDTIDHASHTPGSVS'
               'SAFILEAMVNVISGPKVLMKQIPIWLPLGVADQKTYSFDSTTAAIMLASYTITHFGKATNPLVRVNRLGP'
               'GIPDHPLRLLRIGNQAFLQEFVLPPVQLPQYFTFDLTALKLITQPLPAATWTDDTPTGSNGALRPGISFH'
               'PKLRPILLPNKSGKKGNSADLTSPEKIQAIMTSLQDFKIVPIDPTKNIMGIEVPETLVHKLTGKKVTSKN'
               'GQPIIPVLLPKYIGLDPVAPGDLTMVITQDCDTCHSPASLPAVIEK')

presets_dir = os.path.join(path, 'presets')

def worker(P, recs, kwargs):
    df = P.predict_multiple(recs, **kwargs)
    return df

def get_preset_alleles(name):
    df = pd.read_csv(os.path.join(presets_dir, name+'.csv'),comment='#')
    return list(df.allele)

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

def get_overlapping(index, s, length=9, cutoff=25):
    """Get all mutually overlapping kmers within a cutoff area"""

    g=[s]
    vals = [i for i in range(s, s+cutoff) if i in index]
    for i in range(len(vals)-1):
        if vals[i+1]<=vals[i]+length:
            g.append(vals[i+1])
        else:
            break
    return g

def get_predictor(name='tepitope', **kwargs):
    """Get a predictor"""

    if name == 'netmhciipan':
        return NetMHCIIPanPredictor(**kwargs)
    elif name == 'iedbmhc1':
        return IEDBMHCIPredictor(**kwargs)
    elif name == 'iedbmhc2':
        return IEDBMHCIIPredictor(**kwargs)
    elif name == 'iedbbcell':
        return IEDBBCellPredictor(**kwargs)
    elif name == 'tepitope':
        return TEpitopePredictor(**kwargs)
    elif name == 'mhcflurry':
        return MHCFlurryPredictor(**kwargs)
    else:
        print ('no such predictor %s' %name)
        return

def get_length(data):
    """Get peptide length of a dataframe of predictions"""

    if len(data)>0:
        return len(data.head(1).peptide.max())
    return

'''def get_coords_from_sequence(df, genome, key='peptide'):
    """Get peptide coords from parent protein sequences"""

    def func(x):
        seq = x.translation
        st = seq.find(x[key])
        end = st+len(x[key])
        #print (x['name'],x[key], st, end)
        return pd.Series({'start':st,'end':end})# 'name':x['name'],'peptide':x[key]})

    temp = df.merge(genome[['locus_tag','translation']],
                    left_on='name',right_on='locus_tag')#.set_index(df.index)

    #print (temp[:10])
    temp = temp.apply( lambda r: func(r),1)
    #print (temp[:10])
    df = df.drop(['start','end'],1)
    return df.join(temp)'''

def get_coords(df):
    """Get start end coords from position and length of peptides"""

    if 'start' in df.columns:
        return df
    df['start'] = df.pos.astype(int)
    df['end'] = ( df.pos + df.peptide.str.len() ).astype(int)
    return df

def write_fasta(sequences, id=None, filename='tempseq.fa'):

    if isinstance(sequences, basestring):
        sequences = [sequences]
    out = open(filename, 'w')
    i=1
    for seq in sequences:
        if id == None:
            id='temp%s'%i
        SeqIO.write(SeqRecord(Seq(seq),id,
                    description='temp'), out, 'fasta')
        i+=1
    out.close()
    return filename

def get_sequence(seqfile):
    """Get sequence from fasta file"""

    recs = list(SeqIO.parse(seqfile, 'fasta'))[0]
    sequence = recs.seq.tostring()
    return sequence

def clean_sequence(seq):
    """clean a sequence of invalid characters before prediction"""
    import re
    new = re.sub('[-*_#X]', '', seq)
    return new

def get_nearest(df):
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

def summarize(data):
    """Summarise prediction data"""

    #print 'binders with unique cores: %s' %len(self.getUniqueCores(binders=True))
    allelegrps = data.groupby('allele')
    proteins = data.groupby('name')
    print ('summary: %s peptides in %s proteins and %s alleles' %(len(data),
                                        len(proteins),len(allelegrps)))
    return

def get_cutoffs(pred=None, data=None, cutoff=5):
    """
    Get score cutoffs per allele for supplied predictions at this
    percentile level. Can be used to set global cutoffs for scoring
    arbitary sequences and then check if they are binders.
    """

    q = (1-cutoff/100.) #score quantile value
    cuts={}
    if data is None:
        data = pred.data
    for a,g in data.groupby('allele'):
        cuts[a] = g[pred.scorekey].quantile(q=q)
    cuts = pd.Series(cuts)
    return cuts

def get_standard_mhc1(name):
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

def get_dqp_list(a):
    """Get DRB list in standard format"""

    s = pd.Series(a)
    s = s[s.str.contains('DQ')]
    s = s.apply(lambda x:x.replace('_','*'))
    return list(s)

def get_standard_mhc2(x):
    return 'HLA'+x.replace('_','*')

def compare_predictors(p1, p2, by='allele', cutoff=5, n=2):
    """
    Compare predictions from 2 different predictors.
    Args:
        p1, p2: predictors with prediction results for the same
        set of sequences andalleles
        by: how to group the correlation plots
    """

    import pylab as plt
    import seaborn as sns
    a = p1.promiscuousBinders(n=n, cutoff=cutoff)
    b = p2.promiscuousBinders(n=n, cutoff=cutoff)
    f = utilities.venndiagram([a.peptide, b.peptide],[p1.name,p2.name],colors=('y','b'))
    f.suptitle('common\npromiscuous binders n=%s' %n)
    plt.tight_layout()

    for p in [p1,p2]:
        if not 'score' in p.data.columns:
            p.data['score'] = p.data[p.scorekey]

    #b1 = p1.getBinders(perc=perc)
    #b2 = p2.getBinders(perc=perc)
    #o = analysis.getOverlaps(b1,b2)
    #merge data for plotting score correlation
    df = pd.merge(p1.data, p2.data, on=['peptide','allele','name','pos'])
    g = sns.lmplot(x='score_x', y='score_y', data=df, col=by, ci=None, #robust=True,
               hue=by, col_wrap=3, markers='o', scatter_kws={'alpha':0.5,'s':30})

    g.set_axis_labels(p1.name, p2.name)
    #g.set(xlim=(0, 1))
    #g.set(ylim=(0, 1))
    #plotting.plot_tracks([pi,pf],'MAP_0005',n=2,perc=0.97,legend=True,colormap='jet')
    return

def reshape_data(pred, peptides=None, name=None, values='score'):
    """
    Create summary table per binder/allele with cutoffs applied.
    Args:
        pred: predictor with data
        cutoff: percentile cutoff
        n: number of alleles
    """

    df = pred.data
    idx = ['name','pos','peptide']
    if name != None:
        df = df[df.name==name]
        idx = ['pos','peptide']
    if peptides is not None:
        df = df[df.peptide.isin(peptides)]
    p = df.pivot_table(index=idx, columns='allele', values=values)
    p = p.round(3)
    #get all in p > cutoff for that allele
    order = p.mean(1).sort_values()
    p = p.reindex_axis(order.index, axis=0)
    return p

def plot_summary_heatmap(p, kind='default', name=None):
    """
    Plot heatmap of binders using summary dataframe.
    """

    import pylab as plt
    import seaborn as sns
    sns.set_context("notebook", font_scale=1.2)
    plt.rcParams['savefig.dpi']=100
    h = int(2 + .2 *len(p))
    #print h
    plt.figure(figsize=(10,h))
    ax = sns.heatmap(p, cbar_kws={"shrink": .5})
    #g = sns.clustermap(p, figsize=(10,h), standard_scale=1)#, col_cluster=False)
    #plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    if name != None:
        ax.set_title(name)
    return ax

def protein_summary(pred, peptides, name):
    """formatted protein summary table"""

    x = reshape_data(pred, peptides, name=name)
    return x.style.background_gradient(cmap='Reds', high=.2).highlight_null('white')

def summarize_by_protein(pred, pb):
    """Heatmaps or tables of binders per protein/allele"""

    for i,g in pb.groupby('name'):
        x = binders_allele_summary(pred, g.peptide, values='score', name=i)
        ax=plot_summary_heatmap(x, name=i)

class Predictor(object):
    """Base class to handle generic predictor methods, usually these will
       wrap methods from other modules and/or call command line predictors.
       Subclass for specific functionality"""

    def __init__(self, data=None):
        self.data = data
        self.name = ''
        self.scorekey = 'score'
        self.operator = '<'
        self.rankascending = 0
        #can specify per allele cutoffs here
        self.allelecutoffs = None
        self.temppath = tempfile.mkdtemp()
        return

    def __repr__(self):

        if (self.data is None) or len(self.data) == 0:
            return '%s predictor' %self.name
        else:
            n = len(self.data.name.unique())
            return '%s predictor with results in %s proteins' %(self.name, n)

    def predict(self, sequence, peptide, length=9, overlap=1,
                    allele='', name=''):
        """Does the actual scoring of a sequence. Should be overriden.
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
        return df

    def evaluate(self, df, key, value, operator='<'):
        """
        Evaluate binders less than or greater than a cutoff.
        This method is called by all predictors to get binders
        """

        if operator == '<':
            return df[df[key] <= value]
        else:
            return df[df[key] >= value]

    def getBinders(self, name=None, cutoff=5, cutoff_method='default', data=None, **kwargs):
        """
        Get the top scoring binders. If using scores cutoffs are derived
        from the available prediction data stored in the object. For
        per protein cutoffs the rank can used instead. This will give
        slightly different results.
        Args:
            name: name of protein in predictions, optional
            cutoff: percentile cutoff for score or rank cutoff if value='rank'
            value: 'score' or 'rank'
        Returns:
            binders above cutoff in all alleles, pandas dataframe
        """

        if data is None:
            if self.data is None:
                return
            data = self.data
        cutoff = float(cutoff)
        if name != None:
            if name not in self.proteins():
                print ('no such protein name in binder data')
                return
            data = data[data.name==name]

        if cutoff_method in ['default','']:
            #calculates per allele cutoff over all data laoded
            value = self.scorekey
            if self.rankascending == 0:
                q = (1-cutoff/100.)
            else:
                q = cutoff/100
            if hasattr(self, 'cutoffs'):
                cuts = self.cutoffs
            else:
                cuts={}
                #derive score cutoffs for each allele
                for a,g in data.groupby('allele'):
                    cuts[a] = g[value].quantile(q=q)
            cuts = pd.Series(cuts)
            #print (cuts)
            res=[]
            for a,g in data.groupby('allele'):
                #print cuts[a]
                if self.rankascending == 0:
                    b = g[g[value]>cuts[a]]
                else:
                    b = g[g[value]<cuts[a]]
                res.append(b)
            return pd.concat(res)
        elif cutoff_method == 'rank':
            #done by rank in each sequence/allele
            res = data[data['rank'] < cutoff]
            return res
        elif cutoff_method == 'score':
            #done by global single score cutoff
            #print (data[self.scorekey])
            if self.rankascending == 0:
                res = data[data[self.scorekey] >= cutoff]
            else:
                res = data[data[self.scorekey] <= cutoff]
            return res

    def promiscuousBinders(self, binders=None, name=None, cutoff=5,
                           cutoff_method='default', n=1, unique_core=True, **kwargs):
        """
        Use params for getbinders if no binders provided?
        Args:
            binders: can provide a precalculated list of binders
            name: specific protein, optional
            value: to pass to getBinders
            cutoff: percentile cutoff for getBinders
            n: min number of alleles
            unique_core: removes peptides with duplicate cores and picks the most
            promiscuous and highest ranked, used for mhc-II predictions
        Returns:
            a pandas dataframe
        """

        n=int(n)
        if binders is None:
            binders = self.getBinders(name=name, cutoff=cutoff, cutoff_method=cutoff_method)
        if binders is None:
            print('no binders found, check that prediction data is present')
            return
        if 'core' not in binders.columns :
            binders['core'] = binders.peptide
        grps = binders.groupby(['peptide','pos','name'])
        if self.operator == '<':
            func = min
            skname = 'min'
        else:
            func = max
            skname = 'max'
        s = grps.agg({'allele':pd.Series.count,
                      'core': first, self.scorekey:[func,np.mean],
                      'rank': np.median})
        s.columns = s.columns.get_level_values(1)
        s.rename(columns={skname: self.scorekey, 'count': 'alleles','median':'median_rank',
                         'first':'core'}, inplace=True)
        s = s.reset_index()
        s = s.sort_values(['alleles','median_rank'],ascending=[False,True])

        #if we just want unique cores, drop duplicates takes most promiscuous in each group
        #since we have sorted by alleles and median rank
        if unique_core == True:
            s = s.drop_duplicates('core')
        s = s[s.alleles>=n]
        return s

    def rankedBinders(self, names=None, how='median', cutoff=None):
        """
        Get the median/mean rank of each binder over all alleles.
        Args:
            names: list of protein names, otherwise all current data used
            how: method to use for rank selection, 'median' (default),
            'best' or 'mean',
            cutoff: apply a rank cutoff if we want to filter (optional)
        """

        df = self.data
        if names != None:
            if names is str:
                names = [names]
            df=df[df.name.isin(names)]
        funcs = { 'median':np.median, 'mean':np.mean, 'best':min }
        func = funcs[how]
        b = df.groupby(['peptide']).agg({'rank': func,'pos':first, 'name':first,
                                         self.scorekey: np.median})
        #b.columns = b.columns.get_level_values(0)
        b = b.reset_index().sort_values('rank')
        if cutoff != None:
            b = b[b['rank'] < cutoff]
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

    def predictSequences(self, sequences, alleles=[]):
        """
        Predict a set of arbitary sequences in a list, dict or dataframe.
        These are treated as individual peptides and not split into n-mers.
        """

        results=[]
        #print (sequences)
        if type(sequences) is list:
            sequences = pd.DataFrame(sequences, columns=['peptide'])
            sequences['name'] = sequences.peptide
        for i,row in sequences.iterrows():
            seq = row.peptide
            name = row['name']
            res=[]
            for a in alleles:
               df = self.predict(sequence=seq, length=len(seq),
                                    allele=a, name=name)
               res.append(df)
            res = pd.concat(res)
            results.append(res)
        data = pd.concat(results)
        data.reset_index(drop=True,inplace=True)
        #rank is now just global over all sequences per allele
        data = data.groupby('allele').apply(self.getRanking).reset_index(drop=True)
        self.data = data
        return data

    def predictProteins(self, recs, key='locus_tag', seqkey='translation',
                        names=None, alleles=[], path=None, verbose=False,
                        cpus=1, **kwargs):
        """
        Get predictions for a set of proteins and/or over multiple alleles.
        This is mostly a wrapper for predict_multiple. Sequences should be put into
        a dataframe and supplied to this method as the first argument.
          Args:
            recs: protein sequences in a pandas DataFrame
            key: seq/protein name key
            seqkey: key for sequence column
            names: names of proteins to use from sequences, list or pandas series
            cpus: number of threads to run, use 0 for all cpus
            see predict_multiple for other kwargs
          Returns:
            a dataframe of predictions over multiple proteins
        """

        if type(alleles) is str:
            alleles = [alleles]
        elif type(alleles) is pd.Series:
            alleles = alleles.tolist()
        #self.check_alleles(alleles)
        if len(alleles) == 0:
            print ('no alleles provided')
            return

        if names is not None:
            recs = recs[recs[key].isin(names)]
        results = []
        if path is not None and path != '':
            if not os.path.exists(path):
                os.mkdir(path)

        if verbose == True:
            self.print_heading()
        if cpus == 1:
            results = self.predict_multiple(recs, path, alleles=alleles, seqkey=seqkey,
                                            key=key, verbose=verbose, **kwargs)
        else:
            results = self._multiprocess_predict(recs, path=path, alleles=alleles, seqkey=seqkey,
                                            key=key, verbose=verbose, cpus=cpus, **kwargs )

        print ('predictions done for %s sequences in %s alleles' %(len(recs),len(alleles)))
        if path is None:
            #if no path we keep assign results to the data object
            #assumes we have enough memory..
            self.data = results
        else:
            print ('results saved to %s' %os.path.abspath(path))
            results = None
        self.cleanup()
        return results

    def predict_multiple(self, recs, path=None, overwrite=True, alleles=[], length=11, overlap=1,
                          key='locus_tag', seqkey='sequence', verbose=False,
                          method=None):
        """Predictions for multiple proteins in a dataframe
            Args:
                recs: protein sequences in a pandas DataFrame
                path: if results are to be saved to disk provide a path, otherwise results
                overwrite: over write existing protein files in path if present
                alleles: allele list
                length: length of peptides to predict
                overlap: overlap of n-mers
                key: seq/protein name key
                seqkey: key for sequence column
                verbose: provide output per protein/sequence
                method: IEDB method if using those predictors
            Returns: a dataframe of the results if no path is given
        """

        results = []
        self.length = length
        for i,row in recs.iterrows():
            seq = row[seqkey]
            seq = clean_sequence(seq) #clean the sequence of non-aa characters
            #print (row)
            name = row[key]
            if path is not None:
                fname = os.path.join(path, name+'.csv')
                if os.path.exists(fname) and overwrite == False:
                    continue
            res = []
            for a in alleles:
                #print (a)
                df = self.predict(sequence=seq, length=length, overlap=overlap,
                                    allele=a, name=name, method=method)
                if df is not None:
                    res.append(df)
                else:
                    continue
                if verbose == True and len(df)>0:
                    x = df.iloc[0]
                    s = self.format_row(x)
                    print (s)
            if len(res) == 0:
                continue
            res = pd.concat(res)

            if path is not None and len(res)>0:
                res.to_csv(fname)
            else:
                results.append(res)
        #print (len(recs))
        if len(results)>0:
            results = pd.concat(results)
        return results

    def print_heading(self):
        s = ("{:<30} {:<16} {:<18} {:<}"
                           .format('name','allele','top peptide','score'))
        print (s)
        return

    def format_row(self, x):
        s = ("{:<30} {:<16} {:<18} {:} "
                           .format(x['name'], x.allele, x.peptide, x[self.scorekey] ))
        return s

    def _multiprocess_predict(self, recs, names=[], cpus=2, **kwargs):
        """Call predictproteins with multiprocessing pools
           for running predictions in parallel."""

        import multiprocessing as mp
        maxcpu = mp.cpu_count()
        if cpus == 0 or cpus > maxcpu:
            cpus = maxcpu
        pool = mp.Pool(cpus)
        funclist = []
        st = time.time()
        chunks = np.array_split(recs,cpus)
        for recs in chunks:
            f = pool.apply_async(worker, [self,recs,kwargs])
            #print (f)
            funclist.append(f)
        result = []
        for f in funclist:
            df = f.get(timeout=None)
            if df is not None and len(df)>0:
                result.append(df)
        pool.close()
        pool.join()
        #print (result)
        if len(result)>0:
            result = pd.concat(result)
            #print result.info()
        t=time.time()-st
        print ('took %s' %str(t))
        return result

    def load(self, path=None, names=None,
               compression='infer', file_limit=None):
        """
        Load results for one or more proteins
        Args:
            path: name of a csv file or directory with one or more csv files
            file_limit: limit to load only the this number of proteins
        """

        if os.path.isfile(path):
            self.data = pd.read_csv(path, index_col=0)
        elif os.path.isdir(path):
            if not os.path.exists(path):
                print('no such path %s' %path)
                return
            files = glob.glob(os.path.join(path, '*.csv'))
            if names is not None:
                names = [n+'.csv' for n in names]
                files = [n for n in files if os.path.basename(n) in names]
            if len(files) == 0:
                #print ('no files to load')
                return
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

    def allele_summary(self, cutoff=5):
        """Allele based summary"""

        b = self.getBinders(cutoff=cutoff)
        s = b.groupby('allele').agg({'peptide':np.size,self.scorekey:[np.median,np.mean]})
        s.columns = s.columns.get_level_values(1)
        return s

    def protein_summary(self):
        print ( self.data.groupby('name').agg({'peptide':np.size}) )

    def proteins(self):
        return list(self.data.name.unique())

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
        plot = plotting.plot_tracks([self], name=name, **kwargs)
        return plot

    def getAlleles(self):
        """Get available alleles - override"""
        return []

    def check_alleles(self, alleles):
        a = self.getAlleles()
        found = list((set(a) & set(alleles)))
        if len(found) == 0:
            return
        else:
            return 1

    def cleanup(self):
        """Remove temp files from predictions"""

        if os.path.exists(self.temppath):
            shutil.rmtree(self.temppath)
        return

class NetMHCIIPanPredictor(Predictor):
    """netMHCIIpan predictor"""

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'netmhciipan'
        self.colnames = ['pos','HLA','peptide','Identity','Pos','Core',
                         '1-log50k(aff)','Affinity','Rank']
        self.scorekey = 'Affinity' #'1-log50k(aff)'
        self.cutoff = 500 #.426
        self.operator = '<'
        self.rankascending = 1

    def readResult(self, res):
        """Read raw results from netMHCIIpan output"""

        data=[]
        res = res.split('\n')[19:]
        ignore=['Protein','pos','Number','']
        for r in res:
            if r.startswith('-'): continue
            row = re.split('\s*',r.strip())[:9]
            if len(row)!=9 or row[0] in ignore:
                continue
            data.append(dict(zip(self.colnames,row)))
        return data

    def prepareData(self, df, name):
        """Prepare netmhciipan results as a dataframe"""

        #df = df.convert_objects(convert_numeric=True)
        df = df.apply( lambda x: pd.to_numeric(x, errors='ignore').dropna())
        df['name'] = name
        df.rename(columns={'Core': 'core','HLA':'allele'}, inplace=True)
        df = df.drop(['Pos','Identity','Rank'],1)
        df = df.dropna()
        df['allele'] = df.allele.apply( lambda x: self.convert_allele_name(x) )
        self.getRanking(df)
        self.data = df
        return

    def runSequence(self, seq, length, allele, name='temp', overlap=1):
        """Run netmhciipan for a single sequence"""

        tempfile = os.path.join(self.temppath, name+'.fa')
        seqfile = write_fasta(seq, id=name, filename=tempfile)
        cmd = 'netMHCIIpan -s -length %s -a %s -f %s' %(length, allele, seqfile)
        #print cmd
        temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        #print (temp)
        rows = self.readResult(temp)
        df = pd.DataFrame(rows)
        return df

    def predict(self, sequence=None, peptides=None, length=11, overlap=1,
                    allele='HLA-DRB1*0101', name='',
                    pseudosequence=None, **kwargs):
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
            res = self.runSequence(sequence, length, allele, name, overlap)
        if len(res)==0:
            return res
        self.prepareData(res, name)
        #print (self.data[self.data.columns[:7]][:5])
        return self.data

    def getAlleles(self):
        """Get available alleles"""

        cmd = 'netMHCIIpan -list'
        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except:
            print('netmhciipan not installed?')
            return []
        alleles = temp.split('\n')[34:]
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
        self.methods = ['ann', 'IEDB_recommended', 'comblib_sidney2008',
                        'consensus', 'smm', 'netmhcpan', 'smmpmbec']
        self.scorekey = 'ic50'
        self.cutoff = 500
        self.operator = '<'
        self.rankascending = 1
        self.iedbmethod = 'IEDB_recommended'
        return

    def predict(self, sequence=None, peptides=None, length=11, overlap=1,
                   allele='HLA-A*01:01', name='', method=None):
        """Use IEDB MHCI python module to get predictions.
           Requires that the iedb MHC tools are installed locally"""

        tempfile = os.path.join(self.temppath, name+'.fa')
        seqfile = write_fasta(sequence, id=name, filename=tempfile)
        #print (seqfile)
        path = iedbmhc1path
        if not os.path.exists(path):
            print ('IEDB mhcI tools not found')
            return
        if method == None: method = 'IEDB_recommended'
        self.iedbmethod = method
        cmd = os.path.join(path,'src/predict_binding.py')
        cmd = cmd+' %s %s %s %s' %(method,allele,length,seqfile)
        #print (cmd)
        from subprocess import Popen, PIPE
        try:
            p = Popen(cmd, stdout=PIPE, shell=True)
            temp,error = p.communicate()
            #print (temp)
        except OSError as e:
            print (e)
            return
        df = self.prepareData(temp, name)
        return df

    def prepareData(self, rows, name):
        """Prepare data from results"""

        df = pd.read_csv(io.BytesIO(rows),sep="\t")
        if len(df)==0:
            print (rows) #should print error string from output
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
        if 'method' not in df.columns:
            df['method'] = self.iedbmethod
        if self.iedbmethod in ['IEDB_recommended','consensus']:
            df['ic50'] = df.filter(regex="ic50").mean(1)
        if not 'score' in df.columns:
            df['score'] = df.ic50.apply( lambda x: 1-math.log(x, 50000))
        self.getRanking(df)
        self.data = df
        #print (df[:10])
        return df

    def getAlleles(self):
        """Get available alleles from model_list file and
            convert to standard names"""

        try:
            afile = os.path.join(iedbmhc1path, 'data/MHCI_mhcibinding20130222/consensus/model_list.txt')
            df = pd.read_csv(afile,sep='\t',names=['name','x'])
            alleles = list(df['name'])
            alleles = sorted(list(set([get_standard_mhc1(i) for i in alleles])))
        except:
            alleles = pd.read_csv(os.path.join(datadir, 'iedb_mhc1_alleles.csv')).allele.values
        return alleles

    def getAlleleData(self):
        if not os.path.exists(iedbmhc1path):
            print ('iedb tools not found')
            return
        c = os.path.join(iedbmhc1path,'src/predict_binding.py')
        for m in self.methods:
            cmd = c + ' %s mhc' %m
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
            print (temp)
        return

class IEDBMHCIIPredictor(Predictor):
    """Using IEDB MHC-II method, requires tools to be installed locally"""

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'iedbmhc2'
        self.scorekey = 'score'
        self.cutoff = 4
        self.operator = '<'
        self.rankascending = 0
        self.methods = ['comblib','consensus3','IEDB_recommended',
                        'NetMHCIIpan','nn_align','smm_align','tepitope']
        self.iedbmethod = 'IEDB_recommended'
        return

    def prepareData(self, rows, name):
        """Read data from raw output"""

        if len(rows) == 0:
            return

        #print (rows)
        df = pd.read_csv(io.BytesIO(rows),sep=r'\t',engine='python',index_col=False)
        #print (df.iloc[0])
        extracols = ['Start','End','comblib_percentile','smm_percentile','nn_percentile',
                     'Sturniolo core',' Sturniolo score',' Sturniolo percentile']
        #df = df.drop(extracols,1)
        df.reset_index(inplace=True)
        df.rename(columns={'index':'pos','Sequence': 'peptide','Allele':'allele'},
                           inplace=True)
        df['core'] = df[df.filter(regex="core").columns[0]]
        df['name'] = name

        if self.iedbmethod == 'IEDB_recommended':
            df['score'] = df.percentile_rank
        else:
            if not 'ic50' in df.columns:
                df['ic50'] = df.filter(regex="ic50").mean(1)
            if not 'score' in df.columns:
                df['score'] = df.ic50.apply( lambda x: 1-math.log(x, 50000))

        self.getRanking(df)
        self.data = df
        return df

    def predict(self, sequence=None, peptides=None, length=15, overlap=None,
                   allele='HLA-DRB1*01:01', method='IEDB_recommended', name=''):
        """Use iedb MHCII python module to get predictions.
           Requires that the IEDB MHC-II tools are installed locally
        """

        self.iedbmethod = method
        tempfile = os.path.join(self.temppath, name+'.fa')
        seqfile = write_fasta(sequence, id=name, filename=tempfile)
        path = iedbmhc2path
        if method == None: method = 'IEDB_recommended'
        if not os.path.exists(path):
            print ('iedb mhcII tools not found')
            return
        cmd = os.path.join(path,'mhc_II_binding.py')
        cmd = cmd+' %s %s %s' %(method,allele,seqfile)
        #print (cmd)
        #print (allele)
        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except:
            print ('allele %s not available?' %allele)
            return
        data = self.prepareData(temp, name)
        return data

    def getAlleles(self):
        if not os.path.exists(iedbmhc2path):
            return
        c = os.path.join(iedbmhc2path,'mhc_II_binding.py')
        for m in self.methods:
            cmd = c + ' %s mhc' %m
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
            print (temp)
        return

class TEpitopePredictor(Predictor):
    """Predictor using TepitopePan QM method"""

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'tepitope'
        self.pssms = tepitope.getPSSMs()
        self.cutoff = 2
        self.operator = '>'
        self.rankascending = 0

    def predict(self, sequence=None, peptides=None, length=9, overlap=1,
                    allele='HLA-DRB1*0101', name='',
                    pseudosequence=None, **kwargs):

        self.sequence = sequence
        allele = allele.replace(':','')
        if not allele in self.pssms:
            #print 'computing virtual matrix for %s' %allele
            m = tepitope.createVirtualPSSM(allele)
            if m is None:
                print ('no such allele', allele)
                return pd.DataFrame()
        else:
            m = self.pssms[allele]
        m = m.transpose().to_dict()
        result = tepitope.getScores(m, sequence, peptides, length, overlap=overlap)
        df = self.prepareData(result, name, allele)
        self.data = df
        #print(df[:5])
        return df

    def getAlleles(self):
        return tepitope.getAlleles()

class IEDBBCellPredictor(Predictor):
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
        res=[]
        for i,row in proteins:
            seq = row['translation']
            name = row['locus_tag']
            #print (name)
            df = self.predict(sequence=seq,name=name)
            res.append(df)
            if save == True:
                #fname = os.path.join(path, name+'.mpk')
                #pd.to_msgpack(fname, res)
                fname = os.path.join(path, name+'.csv')
                df.to_csv(fname)
        self.data = res = pd.concat(res)
        return

class MHCFlurryPredictor(Predictor):
    """
    Predictor using MHCFlurry for MHC-I predictions. Requires you to
    install the python package mhcflurry with dependencies.
    see https://github.com/hammerlab/mhcflurry
    """

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'mhcflurry'
        self.cutoff = 500
        self.operator = '<'
        self.scorekey = 'score'
        self.rankascending = 1
        return

    def predict(self, sequence=None, peptides=None, length=11, overlap=1,
                      allele='HLA-A0101', name='', **kwargs):
        """Uses mhcflurry python classes for prediction"""

        self.sequence = sequence
        from mhcflurry import Class1AffinityPredictor
        predictor = Class1AffinityPredictor.load()
        if peptides == None:
            peptides, s = peptutils.create_fragments(seq=sequence,
                                                    length=length, overlap=overlap)
        df = predictor.predict_to_dataframe(peptides=peptides, allele=allele)
        #print (df[:5])
        df = self.prepareData(df, name, allele)
        self.data = df
        return df

    def prepareData(self, df, name, allele):
        """Post process dataframe to alter some column names"""

        df = df.rename(columns={'Allele':'allele','Peptide':'peptide'})
        df['name'] = name
        df['pos'] = df.index
        #df['score'] = df['prediction'].apply( lambda x: 1-math.log(x, 50000) )
        df['score'] = df.prediction
        self.getRanking(df)
        return df

    def convert_allele_name(self, r):
        return r[:5]+'*'+r[5:7]+':'+r[7:]

    def getAlleles(self):
        import mhcflurry
        return mhcflurry.Class1AffinityPredictor.supported_alleles
