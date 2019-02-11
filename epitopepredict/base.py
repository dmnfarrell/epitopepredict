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
from . import utilities, peptutils, sequtils, tepitope, config, mhclearn

#write a default config if not present
config_file = config.write_default_config()
home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'mhcdata')
predictors = ['basicmhc1','tepitope','netmhciipan','netmhcpan','mhcflurry','mhcnuggets','iedbmhc1','iedbmhc2']
iedbmhc1_methods = ['ann', 'IEDB_recommended', 'comblib_sidney2008', 'consensus', 'smm', 'smmpmbec']
mhc1_predictors = ['basicmhc1','netmhcpan','iedbmhc1','mhcflurry','mhcnuggets'] + iedbmhc1_methods
iedbsettings = {'cutoff_type': 'none', 'pred_method': 'IEDB_recommended',
            'output_format': 'ascii', 'sort_output': 'position_in_sequence',
            'sequence_format': 'auto', 'allele': 'HLA-DRB1*01:01', 'length':'11',
            'sequence_file': None}

mhc1_presets = ['mhc1_supertypes','us_caucasion_mhc1','us_african_mhc1','broad_coverage_mhc1']
mhc2_presets = ['mhc2_supertypes','human_common_mhc2','bovine_like_mhc2']

#sequence for testing
testsequence = ('MRRVILPTAPPEYMEAIYPVRSNSTIARGGNSNTGFLTPESVNGDTPSNPLRPIADDTIDHASHTPGSVS'
               'SAFILEAMVNVISGPKVLMKQIPIWLPLGVADQKTYSFDSTTAAIMLASYTITHFGKATNPLVRVNRLGP'
               'GIPDHPLRLLRIGNQAFLQEFVLPPVQLPQYFTFDLTALKLITQPLPAATWTDDTPTGSNGALRPGISFH'
               'PKLRPILLPNKSGKKGNSADLTSPEKIQAIMTSLQDFKIVPIDPTKNIMGIEVPETLVHKLTGKKVTSKN'
               'GQPIIPVLLPKYIGLDPVAPGDLTMVITQDCDTCHSPASLPAVIEK')
presets_dir = os.path.join(module_path, 'presets')

def read_defaults():
    """Get some global settings such as program paths from config file"""

    cp = config.parse_config(config_file)
    options = config.get_options(cp)
    return options

#get defaults
defaults = read_defaults()

def predict_proteins_worker(P, recs, kwargs):
    try:
        df = P._predict_sequences(recs, **kwargs)
    except KeyboardInterrupt:
        return
    return df

def predict_peptides_worker(P, recs, kwargs):
    try:
        df = P._predict_peptides(recs, **kwargs)
    except KeyboardInterrupt:
        return
    return df

def get_preset_alleles(name):
    df = pd.read_csv(os.path.join(presets_dir, name+'.csv'),comment='#')
    return list(df.allele)

def first(x):
    return x.iloc[0]

def get_pos(x):
    x['pos'] = np.arange(len(x))
    return x

def get_iedb_request(seq, alleles='HLA-DRB1*01:01', method='consensus3'):
    import requests
    url = 'http://tools.iedb.org/tools_api/mhcii/'
    values = {'method' : method,
              'sequence_text' : seq,
              'allele' : alleles }
    r=requests.post(url,data=values)
    df=pd.read_csv(io.StringIO(r.content),sep='\t')
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

def get_predictor_classes():
    """Get predictor classes in this module."""

    cl={}
    for o in Predictor.__subclasses__():
        cl[o._get_name()] = o
    return cl

def get_predictor(name='tepitope', **kwargs):
    """Get a predictor object using it's name. Valid predictor names
       are held in the predictors attribute."""

    cl = get_predictor_classes()
    if name in iedbmhc1_methods:
        name = 'iedbmhc1'
    if name in cl:
        return cl[name](**kwargs)
    else:
        print ('no such predictor %s' %name)
        print ('valid names are %s' %', '.join(predictors))
        return

def get_length(data):
    """Get peptide length of a dataframe of predictions"""

    if len(data)>0:
        return len(data.head(1).peptide.max())
    return

def get_coords(df):
    """Get start end coords from position and length of peptides"""

    if 'start' in df.columns:
        return df
    df['start'] = df.pos.astype(int)
    df['end'] = ( df.pos + df.peptide.str.len() ).astype(int)
    return df

def seq_from_binders(df):
    x = pd.Series(df.sort_values('pos').peptide.unique())
    return ''.join(x.str[0])

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

def get_filenames(path, names=None, file_limit=None):

    if not os.path.exists(path):
        print('no such path %s' %path)
        return
    files = glob.glob(os.path.join(path, '*.csv'))
    if names is not None:
        names = [n+'.csv' for n in names]
        files = [n for n in files if os.path.basename(n) in names]
    if len(files) == 0:
        return
    if file_limit != None:
        files = files[:file_limit]
    return files

def results_from_csv(path=None, names=None, compression='infer', file_limit=None):
    """
    Load results for multiple csv files in a folder or a single file.
    Args:
        path: name of a csv file or directory with one or more csv files
        names: names of proteins to load
        file_limit: limit to load only the this number of proteins
    """

    if os.path.isfile(path):
        data = pd.read_csv(path, index_col=0)
    elif os.path.isdir(path):
        files = get_filenames(path, names, file_limit)
        if files == None:
            return
        res = []
        for f in files:
            df = pd.read_csv(f, index_col=0, compression=compression)
            if len(df) == 0:
                continue
            #if not self.scorekey in df.columns:
            #    continue
            res.append(df)
        data = pd.concat(res)
    else:
        return
    return data

def get_quantiles(predictor):
    """Get quantile score values per allele in set of predictions.
    Used for making pre-defined cutoffs.
    Args:
        predictor: predictor with set of predictions
    """

    P=predictor
    res=[]
    value = P.scorekey
    data=P.data
    for a,g in data.groupby('allele'):
        q = g[value].quantile(np.arange(0,1,.01))
        q.name=a
        res.append(q)
    res = pd.DataFrame(res)
    if P.rankascending==1:
        res.columns = res.columns.sort_values(ascending=False)
    return res

def get_standard_mhc1(name):
    """Taken from iedb mhc1 utils.py"""

    temp = name.strip().split('-')
    length = temp[-1]
    mhc = '-'.join(temp[0:-1])
    return mhc

def get_drb_list(a):
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
    a = p1.promiscuous_binders(n=n, cutoff=cutoff)
    b = p2.promiscuous_binders(n=n, cutoff=cutoff)
    f = utilities.venndiagram([a.peptide, b.peptide],[p1.name,p2.name],colors=('y','b'))
    f.suptitle('common\npromiscuous binders n=%s' %n)
    plt.tight_layout()

    for p in [p1,p2]:
        if not 'score' in p.data.columns:
            p.data['score'] = p.data[p.scorekey]

    #b1 = p1.get_binders(perc=perc)
    #b2 = p2.get_binders(perc=perc)
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

def split_peptides(df,length=9,seqkey='sequence',newcol='peptide'):
    """Split sequences in a dataframe into peptide fragments"""

    res=[]
    for n,r in df.iterrows():
        seq = r[seqkey]
        p = peptutils.get_fragments(seq,1,length)
        p = p.rename(columns={'peptide':newcol})
        p['name']=n
        p[seqkey] = seq
        p=p.set_index('name')
        res.append(p)
    res = pd.concat(res)
    res = df.merge(res,on=seqkey)
    return res

def check_snap():
    """Check if inside a snap"""

    if 'SNAP_COMMON' in os.environ:
        return True
    return False

def set_netmhcpan_cmd(path=None):
    """Setup the netmhcpan command to point directly to the binary. This is a
    workaround for running inside snaps. Avoids using the tcsh script."""

    toolspath = os.path.join('/home', os.environ['USER'], 'tools')
    netmhcpath = os.path.join(toolspath, 'netMHCpan-4.0/Linux_x86_64')
    if not os.path.exists(netmhcpath):
        print ('you need to copy the netmhcpan files to %s' %toolspath)
    os.environ['NETMHCpan']=netmhcpath
    os.environ['TMPDIR']= '/tmp'
    cmd = os.path.join(netmhcpath, 'bin/netMHCpan')
    #print (os.environ)
    return cmd

class DataFrameIterator:
    """Simple iterator to get dataframes from a path out of memory"""
    def __init__(self, files):
        self.files = files
        self.current = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.current >= len(self.files):
            raise StopIteration
        else:
            df=pd.read_csv(self.files[self.current],index_col=0)
            self.current += 1
            return df

    def next(self):
        return self.__next__()

    def __repr__(self):
        return 'DataFrameIterator with %s files' %len(self.files)

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
        self.path = None
        return

    @classmethod
    def _get_name(self):
        return 'base class'

    def __repr__(self):

        if (self.data is None) or len(self.data) == 0:
            return '%s predictor' %self.name
        else:
            n = len(self.data.name.unique())
            return '%s predictor with results in %s sequences' %(self.name, n)

    def check_install(self):
        return True

    def supported_lengths(self):
        """Return supported peptide lengths"""
        return [9,11]

    def predict(self, sequence=None, peptides=None, length=9, overlap=1,
                    allele='', name=''):
        """Does the actual scoring of a sequence. Should be overriden.
           Should return a pandas DataFrame"""
        return

    def prepare_data(self, result, name, allele):
        """Put raw prediction data into DataFrame and rank,
           override for custom processing. Can be overriden for
           custom data."""

        df = pd.DataFrame(result, columns=['peptide','core','pos','score'])
        df['name'] = name
        df['allele'] = allele
        self.get_ranking(df)
        return df

    def get_ranking(self, df):
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

    def get_quantile_data(self):
        """Get peptide score rank from quantile data."""

        path = os.path.join(datadir, 'quantiles_%s.csv' %self.name)
        qf = pd.read_csv(path,index_col=0)
        qf.columns = qf.columns.astype('float')
        return qf

    def get_global_rank(self, score, allele):
        """Get an allele specific score percentile from precalculated quantile data."""

        qf = self.qf
        s = qf.loc[allele]
        cutoff = np.abs(s - score).idxmin()
        return cutoff

    def get_allele_cutoffs(self, cutoff=.95):
        """Get per allele percentile cutoffs using precalculated quantile vales."""

        qf = self.qf
        cutoff = round(cutoff,2)
        #if self.rankascending == 1:
        #    cutoff = round(1-cutoff,2)
        if not cutoff in qf.columns:
            print ('please use values between 0 and 1 rounded to 2 decimal places e.g. .95 (95% cutoff)')
        return qf[cutoff]

    def _per_allele_binders(self, data, cuts):
        """Return binders per allele based cutoffs"""

        value = self.scorekey
        res=[]
        for a,g in data.groupby('allele'):
            if a not in cuts:
                print ('%s not in cutoffs' %a)
                continue
            #print (cuts[a])
            #print (g[value])
            if self.rankascending == 0:
                b = g[g[value]>=cuts[a]]
            else:
                b = g[g[value]<cuts[a]]
            res.append(b)
        return pd.concat(res)

    def _ranked_binders(self, data, cutoff):
        """return binders by rank which has been calculated per sequence/allele"""
        return data[data['rank'] <= cutoff]

    def _score_binders(self, data, cutoff):
        """return binders by global single score cutoff"""
        if self.rankascending == 0:
            res = data[data[self.scorekey] >= cutoff]
        else:
            res = data[data[self.scorekey] <= cutoff]
        return res

    def get_scores(self, allele):
        """Return peptides and scores only for an allele"""

        if self.data is None or len(self.data)==0:
            return
        df = self.data[self.data.allele==allele]
        sc = df[['peptide','score']]
        return sc

    def get_binders(self, cutoff=.95, cutoff_method='default', path=None,
                    name=None, drop_columns=False, **kwargs):
        """
        Get the top scoring binders. If using default cutoffs are derived
        from the pre-defined percentile cutoffs for some known antigens. For
        per protein cutoffs the rank can used instead. This will give
        slightly different results.
        Args:
            path: use results in a path instead of loading at once, conserves memory
            cutoff: percentile cutoff (default), absolute score or a rank value within each sequence
            cutoff_method: 'default', 'score' or 'rank'
            name: name of a specific protein/sequence
        Returns:
            binders above cutoff in all alleles, pandas dataframe
        """

        if cutoff_method not in ['default', 'global', 'score', 'rank']:
            print ('choose a valid cutoff method')
            return
        cutoff = float(cutoff)
        data = self.data
        if cutoff_method in ['default','global','']:
            #by per allele percentile cutoffs
            cuts = self.get_allele_cutoffs(cutoff)
        #print (cuts)
        if path is not None:
            #get binders out of memory for large datasets
            files = get_filenames(path)
            if files == None:
                return

            if name is not None:
                files = [f for f in files if f.find(name+'.csv')!=-1]
            d = DataFrameIterator(files)
            res=[]
            for data in d:
                if cutoff_method in ['default','global','']:
                    b = self._per_allele_binders(data, cuts)
                elif cutoff_method == 'rank':
                    b = self._ranked_binders(data, cutoff)
                elif cutoff_method == 'score':
                    b = self._score_binders(data, cutoff)
                res.append(b)
            res = pd.concat(res)
        else:
            if data is None:
                return
            if cutoff_method in ['default','global','']:
                res = self._per_allele_binders(data, cuts)
            elif cutoff_method == 'rank':
                res = self._ranked_binders(data, cutoff)
            elif cutoff_method == 'score':
                res = self._score_binders(data, cutoff)
            names = list(res['name'].unique())
            if name != None and name in names:
                res = res[res.name == name]

        return res

    def promiscuous_binders(self, binders=None, name=None, cutoff=.95,
                           cutoff_method='default', n=1, unique_core=True, **kwargs):
        """
        Use params for getbinders if no binders provided?
        Args:
            binders: can provide a precalculated list of binders
            name: specific protein, optional
            value: to pass to get_binders
            cutoff_method: 'default', 'score' or 'rank'
            cutoff: percentile cutoff for get_binders
            n: min number of alleles
            unique_core: removes peptides with duplicate cores and picks the most
            promiscuous and highest ranked, used for mhc-II predictions
        Returns:
            a pandas dataframe
        """

        n=int(n)
        if binders is None:
            binders = self.get_binders(name=name, cutoff=cutoff, cutoff_method=cutoff_method)

        if binders is None:
            print('no binders found, check that prediction data is present')
            return
        if 'core' not in binders.columns :
            binders['core'] = binders.peptide
        if self.operator == '<':
            func = min
            skname = 'min'
        else:
            func = max
            skname = 'max'

        s = binders.groupby(['peptide','pos','name']).agg({'allele': pd.Series.count,
                            'core': first, self.scorekey:[func,np.mean],
                            'rank': np.median})
        s.columns = s.columns.get_level_values(1)
        s.rename(columns={skname: self.scorekey, 'count': 'alleles','median':'median_rank',
                         'first':'core'}, inplace=True)
        s = s.reset_index()

        '''if keep_columns == True:
            x = binders.drop_duplicates(['name','peptide','allele'])
            s = s.merge(x, on=['name','peptide','pos'], how='inner', suffixes=('', '_y'))
            cols = s.columns[~s.columns.str.contains('_y')]
            s = s[cols]
            s = s.drop('allele',1)'''

        s = s.sort_values(['alleles','median_rank',self.scorekey],
                          ascending=[False,True,self.rankascending])
        #if we just want unique cores, drop duplicates takes most promiscuous in each group
        #since we have sorted by alleles and median rank
        if unique_core == True:
            s = s.drop_duplicates('core')
        s = s[s.alleles>=n]
        return s

    def ranked_binders(self, names=None, how='median', cutoff=None):
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

    def get_unique_cores(self, binders=False):
        """Get only unique cores"""

        if binders == True:
            df = self.get_binders()
        else:
            df = self.data
        grouped = df.groupby('core')
        cores = grouped.agg({self.scorekey:max})
        #cores = df.loc[grouped[self.scorekey].max().index]
        cores.sort(self.scorekey, inplace=True, ascending=self.rankascending)
        #print cores
        return cores

    def seqs_to_dataframe(self, seqs):
        df = pd.DataFrame(seqs, columns=['peptide'])
        df['name'] = df.peptide
        df['pos'] = range(len(seqs))
        return df

    def _predict_peptides(self, peptides, alleles=[], verbose=False, drop_columns=False, **kwargs):
        """
        Predict a set of arbitary peptide sequences in a list or dataframe.
        These are treated as individual peptides and not split into n-mers. This is
        usually called wrapped by predict_peptides. If the predict method for the class
        can only accept protein sequences this needs to be overriden.
        Args:
            peptides: list of peptides
            alleles: list of alleles to predict
            drop_columns: only keep default columns
        Returns:
            dataframe with results
        """

        res=[]
        if len(alleles)==0:
            print ('no alleles specified')
            return
        if type(alleles) is str:
            alleles = [alleles]

        #remove empty values
        #peptides = [i for i in peptides if type(i) is str or type(i) is unicode]
        peptides = [i for i in peptides if isinstance(i, basestring)]
        #print (peptides)

        for a in alleles:
            df = self.predict(peptides=peptides, allele=a, **kwargs)
            if df is None or len(df) == 0:
                continue
            #if keep_empty == True:
                 #df = s.merge(df,on=['peptide'],how='left')
                 #df['allele'] = df.allele.ffill()
                 #df['pos'] = df.index.copy()
            res.append(df)
            if verbose == True and len(df)>0:
                x = df.iloc[0]
                r = self.format_row(x)
                print (r)

        if len(res) == 0:
            return

        data = pd.concat(res)
        #drop non-standard columns for cleaner output
        if drop_columns == True:
            defaultcols = ['pos','peptide',self.scorekey,'allele','rank','name']
            data= data[defaultcols]
        #print (data)
        self.data = data
        return data

    def predict_peptides(self, peptides, cpus=1, path=None, overwrite=True, name=None, **kwargs):
        """Predict a set of individual peptides without splitting them.
        This is a wrapper for _predict_peptides to allow multiprocessing.
        Args:
            peptides: list of peptides
            alleles: list of alleles to predict
            drop_columns: only keep default columns
        Returns:
            dataframe with results
        """

        if path is not None:
            fname = os.path.join(path, 'results_%s.csv' %self.name)
            if overwrite == False and os.path.exists(fname):
                #load the data if present and not overwriting
                print ('results file found, not overwriting')
                self.path = fname
                return

        #store original peptide list as dataframe to keep order
        s = pd.DataFrame(peptides, columns=['peptide'])
        s['pos'] = s.index.copy()

        if cpus == 1:
            data = self._predict_peptides(peptides, **kwargs)
        else:
            data = self._run_multiprocess(peptides, worker=predict_peptides_worker,
                                          cpus=cpus, **kwargs)
        if data is None or len(data) == 0:
            print ('empty result returned')
            return

        #get original peptide positions by re-merging with original peptide dataframe
        #also allows us to keep empty rows in correct order
        #make this optional as it could cause issues?
        new = []
        for i,g in data.groupby('allele'):
            #print (len(g))
            #print (s)
            #print (g)
            if 'pos' in g.columns:
                g=g.drop('pos',1)
            x = s.merge(g,on=['peptide'],how='left').drop_duplicates(['pos','peptide'])
            #x['pos'] = x.index.copy()
            x['allele'] = x.allele.fillna(i)
            new.append(x)
            #print (len(x))
            #print (x)

        data = pd.concat(new).reset_index(drop=True)
        data = data.groupby('allele').apply(self.get_ranking) #move this to loop
        data = data.reset_index(drop=True)

        if name==None:
            name='temp'
        if 'name' not in data.columns:
            data['name'] = name
        #print (data)
        if path is not None:
            data.to_csv(fname, float_format='%g')
            self.path = fname
        else:
            self.data = data
        return data

    def _convert_to_dataframe(self, recs):
        """convert sequence list input to dataframe"""
        if type(recs) is str:
            recs = [recs]
        if type(recs) is list:
            idx = [str(i) for i in range(len(recs))]
            recs = pd.DataFrame(zip(idx,recs), columns=['locus_tag','translation'])
        return recs

    def predict_proteins(self, args, **kwargs):
        """Alias to predict_sequences"""
        results = self.predict_sequences(args, **kwargs)
        return results

    def predict_sequences(self, recs, alleles=[], path=None, verbose=False,
                          names=None, key='locus_tag', seqkey='translation', cpus=1, **kwargs):
        """
        Get predictions for a set of proteins over multiple alleles that allows
        running in parallel using the cpus parameter.
        This is a wrapper for _predictSequences with the same args.
          Args:
            recs: list or dataframe with sequences
            cpus: number of processors
            key: seq/protein name key
            seqkey: key for sequence column
          Returns:
            a dataframe of predictions over multiple proteins
        """

        recs = self._convert_to_dataframe(recs)
        if names is not None:
            recs = recs[recs[key].isin(names)]
        if verbose == True:
            self.print_heading()
        if cpus == 1:
            results = self._predict_sequences(recs, alleles=alleles, path=path, verbose=verbose, **kwargs)
        else:
            results = self._run_multiprocess(recs, alleles=alleles, path=path, verbose=verbose, cpus=cpus, **kwargs)

        print ('predictions done for %s sequences in %s alleles' %(len(recs),len(alleles)))
        if path is None:
            #if no path we assign results to the data attribute
            self.data = results
        else:
            print ('results saved to %s' %os.path.abspath(path))
            results = None
        self.cleanup()
        return results

    def _predict_sequences(self, recs, path=None, overwrite=True, alleles=[], length=11, overlap=1,
                          verbose=False, compression=None, **kwargs):
        """
        Get predictions for a set of proteins and/or over multiple alleles.
        Sequences should be put into
        a dataframe and supplied to this method as the first argument.
            Args:
                recs: protein sequences in a pandas DataFrame
                path: if results are to be saved to disk provide a path, otherwise results
                overwrite: over write existing protein files in path if present
                alleles: allele list
                length: length of peptides to predict
                overlap: overlap of n-mers
                verbose: provide output per protein/sequence
            Returns: a dataframe of the results if no path is given
        """

        if type(alleles) is str:
            alleles = [alleles]
        elif type(alleles) is pd.Series:
            alleles = alleles.tolist()
        #self.check_alleles(alleles)
        if len(alleles) == 0:
            print ('no alleles provided')
            return

        #recs = self._convert_to_dataframe(recs)
        results = []
        if path is not None and path != '':
            if not os.path.exists(path):
                try:
                    os.mkdir(path)
                except (OSError, FileExistsError):
                    pass

        results = []
        self.length = length
        if compression == '':
            compression = None
        for i,row in recs.iterrows():
            seq = row.translation
            seq = clean_sequence(seq) #clean the sequence of non-aa characters
            #print (row)
            name = row.locus_tag
            if path is not None:
                fname = os.path.join(path, name+'.csv')
                if os.path.exists(fname) and overwrite == False:
                    continue
            res = []
            peptides, s = peptutils.create_fragments(seq=seq, length=length, overlap=overlap)
            for a in alleles:
                #we pass sequence, length, overlap for predictors that have to run on whole seq
                df = self.predict(peptides=peptides, sequence=seq, allele=a, name=name,
                                  length=length, overlap=overlap, **kwargs)
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
            res = pd.concat(res, sort=True)

            if path is not None and len(res)>0:
                if compression == 'gzip':
                    fname = fname+'.gz'
                res.to_csv(fname, compression=compression)
            else:
                results.append(res)
        #print (len(recs))
        if len(results)>0:
            results = pd.concat(results)
        return results

    def print_heading(self):
        s = ("{:<25} {:<16} {:<18} {:<}"
                           .format('name','allele','top peptide','score'))
        print (s)
        return

    def format_row(self, x):
        s = ("{:<25} {:<16} {:<18} {:} "
                           .format(x['name'], x.allele, x.peptide, x[self.scorekey] ))
        return s

    def _run_multiprocess(self, recs, cpus=2, worker=None, **kwargs):
        """
        Call function with multiprocessing pools. Used for running predictions
        in parallel where the main input is a pandas dataframe.
        Args:
            recs: input dataframe
            cpus: number of cores to use
            worker: function to be run in parallel
        Returns:
            concatenated result, a pandas dataframe
           """

        if worker is None:
            worker = predict_proteins_worker
            #print (worker)
        import multiprocessing as mp
        maxcpu = mp.cpu_count()
        if cpus == 0 or cpus > maxcpu:
            cpus = maxcpu
        if cpus >= len(recs):
            cpus = len(recs)
        pool = mp.Pool(cpus)
        funclist = []
        st = time.time()
        chunks = np.array_split(recs,cpus)
        #print ([len(i) for i in chunks])

        for recs in chunks:
            f = pool.apply_async(worker, [self,recs,kwargs])
            funclist.append(f)
        result = []
        try:
            for f in funclist:
                df = f.get(timeout=None)
                #print (df)
                if df is not None and len(df)>0:
                    result.append(df)
        except KeyboardInterrupt:
            print ('process interrupted')
            pool.terminate()
            sys.exit(0)
        pool.close()
        pool.join()

        #print (result)
        if len(result)>0:
            result = pd.concat(result)
            #print result.info()
        t=time.time()-st
        if t>10:
            print ('took %s seconds' %str(round(t,3)))
        return result

    def load(self, path=None, names=None, compression='infer', file_limit=None):
        """Load results from path or single file. See results_from_csv for args."""

        if path is None and self.path != None:
            #if path attribute is set no path arg we use that
            path = self.path
        data = results_from_csv(path, names, compression, file_limit)
        if data is None:
            print ('no data found in path %s' %path)
        else:
            self.data = data
        if not self.scorekey in data.columns:
            print ('WARNING. this does not look like the correct data for this predictor')
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

        b = self.get_binders(cutoff=cutoff)
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

    def get_names(self):
        """Get names of sequences currently stored as predictions"""

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

    def get_alleles(self):
        """Get available alleles - override"""
        return []

    def check_alleles(self, alleles):
        a = self.get_alleles()
        #print (a)
        found = list((set(a) & set(alleles)))
        return found

    def cleanup(self):
        """Remove temp files from predictions"""

        if os.path.exists(self.temppath):
            shutil.rmtree(self.temppath)
        return

    def check_install(self):
        return True

class NetMHCPanPredictor(Predictor):
    """netMHCpan 4.0 predictor
    see http://www.cbs.dtu.dk/services/NetMHCpan/
    Default scoring is ligand likelihood predictions.
    To get newer scoring behaviour pass scoring='ligand' to constructor.
    """
    def __init__(self, data=None, scoring='affinity'):

        Predictor.__init__(self, data=data)
        self.name = 'netmhcpan'
        self.colnames = ['pos','allele','peptide','core','Of','Gp','Gl','Ip','Il','Icore',
                         'name','raw_score','ic50','rank']
        self.scorekey = 'score'
        self.scoring = scoring
        if self.scoring =='affinity':
            self.cutoff = 500
            self.operator = '<'
            self.rankascending = 1
        elif scoring == 'ligand':
            self.cutoff = 0.5
            self.operator = '>'
            self.rankascending = 0
        #load precalculated quantiles for sample peptides
        self.qf = self.get_quantile_data()
        self.basecmd = 'netMHCpan'
        #base command needs to be to the binary directly if running snap
        if check_snap() is True:
            self.basecmd = set_netmhcpan_cmd()

    def read_result(self, temp):
        """Read raw results from netMHCpan output"""

        ignore = ['Pos','#','Protein','']
        if self.scoring == 'affinity':
            cols = ['pos','allele','peptide','core','Of','Gp','Gl','Ip','Il','Icore',
                    'Identity','Score','ic50','%Rank','X','BindLevel']
        else:
            cols = ['pos','allele','peptide','core','Of','Gp','Gl','Ip','Il','Icore',
                    'Identity','Score','%Rank','X','BindLevel']
        res = pd.read_csv(io.BytesIO(temp), comment='#', names=cols, sep='\s+',
                          error_bad_lines=False, skiprows=45).dropna(subset=['peptide'])
        res = res.drop(columns='X')
        res = res[~res.pos.isin(ignore)]
        return res

    def prepare_data(self, df, name):
        """Prepare netmhcpan results"""

        df = df.apply( lambda x: pd.to_numeric(x, errors='ignore').dropna())
        df['name'] = name
        if self.scoring =='affinity':
            df['score'] = df.ic50
        else:
            df['score'] = df.Score
        df=df.reset_index(drop=True)
        df['pos'] = df.index.copy()
        df = df.dropna(how='all')
        self.get_ranking(df)
        self.data = df
        return

    def predict(self, peptides, allele='HLA-A*01:01', name='temp',
                pseudosequence=None, show_cmd=False, **kwargs):
        """Call netMHCpan command line.  """

        allele = allele.replace('*','')
        #write peptides to temp file
        pepfile = tempfile.mktemp()+'.pep'
        with open(pepfile ,'w') as f:
            for p in peptides:
                f.write(p+'\n')
        f.close()
        if self.scoring =='affinity':
            cmd = '%s -BA -f %s -inptype 1 -a %s' %(self.basecmd, pepfile , allele)
        else:
            cmd = '%s -f %s -inptype 1 -a %s' %(self.basecmd, pepfile , allele)
        if show_cmd is True:
            print (cmd)
        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except Exception as e:
            #print (e)
            return
        res = self.read_result(temp)

        if len(res)==0:
            print ('empty result, check allele %s is correct' %allele)
            return res
        self.prepare_data(res, name=name)
        return self.data

    def get_alleles(self):
        """Get available alleles"""

        cmd = 'netMHCpan -listMHC'
        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except:
            print('netmhcpan not installed?')
            return []
        s = temp.decode().split('\n')
        alleles = [i for i in s if not i.startswith('#') and len(i)>0]
        alleles = [self.convert_allele_name(i) for i in alleles]
        return alleles

    def convert_allele_name(self, a):
        """Convert allele names to internally used form"""

        if not ':' in a and a.startswith('HLA'):
            a = a[:-2]+':'+a[-2:]
        return a

    def check_install(self):

        cmd = '%s -h' %self.basecmd
        try:
            temp = subprocess.check_output(cmd, shell=True)
            print('netMHCpan appears to be installed')
            return 1
        except:
            return 0

    @classmethod
    def _get_name(self):
        return 'netmhcpan'

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
        #load precalculated quantiles for sample peptides
        self.qf = self.get_quantile_data()

    def read_result(self, res):
        """Read raw results from netMHCIIpan output"""

        data=[]
        res = res.decode()
        res = res.split('\n')
        ignore=['Protein','pos','Number','#']
        for r in res:
            if r.startswith('-'): continue
            row = re.split('\s*',r.strip())[:9]
            if len(row)!=9 or row[0] in ignore:
                continue
            #print (row)
            data.append(dict(zip(self.colnames,row)))
        return data

    def prepare_data(self, df, name):
        """Prepare netmhciipan results as a dataframe"""

        #df = df.convert_objects(convert_numeric=True)
        df = df.apply( lambda x: pd.to_numeric(x, errors='ignore').dropna())
        df['name'] = name
        df.rename(columns={'Core': 'core','HLA':'allele'}, inplace=True)
        df = df.drop(['Pos','Identity','Rank'],1)
        df['allele'] = df.allele.apply( lambda x: self.convert_allele_name(x) )
        df['score'] = pd.to_numeric(df.Affinity,errors='coerce')
        df = df.dropna()
        self.get_ranking(df)
        self.data = df
        return

    def predict(self, peptides, allele='HLA-DRB1*0101', name='temp',
                pseudosequence=None, show_cmd=False, **kwargs):
        """Call netMHCIIpan command line."""

        #assume allele names are in standard format HLA-DRB1*0101
        #if allele.startswith('HLA-'):
        #    allele=allele[4:]
        allele = self.allele_mapping(allele)
        #print (allele)
        #write peptides to temp file
        pepfile = tempfile.mktemp()+'.pep'
        with open(pepfile ,'w') as f:
            for p in peptides:
                f.write(p+'\n')
        f.close()

        cmd = 'netMHCIIpan -f %s -inptype 1 -a %s' %(pepfile , allele)
        if show_cmd is True:
            print (cmd)
        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except Exception as e:
            #print (e)
            return
        rows = self.read_result(temp)
        res = pd.DataFrame(rows)
        if len(res)==0:
            return res
        self.prepare_data(res, name)
        return self.data

    def allele_mapping(self, allele):
        if allele.startswith('HLA-DQA'):
            return allele.replace('*','').replace(':','')
        else:
            return allele.replace('*','_').replace(':','').replace('HLA-DRB','DRB')

    def get_alleles(self):
        """Get available alleles"""

        cmd = 'netMHCIIpan -list'
        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except:
            print('netmhciipan not installed?')
            return []
        alleles = temp.decode().split('\n')#[34:]
        return alleles

    def convert_allele_name(self, a):
        """Convert allele names to internally used form"""

        if not a.startswith('HLA'):
            return 'HLA-'+a.replace(':','')
        else:
            return a.replace(':','')

    def check_install(self):
        cmd = 'netMHCIIpan'
        try:
            temp = subprocess.check_output(cmd, shell=True)
            print('netMHCIIpan appears to be installed')
            return 1
        except:
            print ('netMHCIIpan not found')
            return 0

    @classmethod
    def _get_name(self):
        return 'netmhciipan'

class IEDBMHCIPredictor(Predictor):
    """Using IEDB tools method, requires iedb-mhc1 tools. Tested with version 2.17"""

    def __init__(self, data=None, method='IEDB_recommended'):
        Predictor.__init__(self, data=data)
        self.name = 'iedbmhc1'
        self.methods = ['ann', 'IEDB_recommended', 'comblib_sidney2008',
                        'consensus', 'smm', 'netmhcpan', 'smmpmbec']
        self.scorekey = 'ic50'
        self.cutoff = 500
        self.operator = '<'
        self.rankascending = 1
        self.method = method
        self.iedbmhc1_path = defaults['iedbmhc1_path']
        self.qf = self.get_quantile_data()
        return

    def predict(self, sequence=None, peptides=None, length=11, overlap=1,
                   allele='HLA-A*01:01', name='', method=None, show_cmd=False, **kwargs):
        """Use IEDB MHCI python module to get predictions.
           Requires that the IEDB MHC tools are installed locally
           Args:
            sequence: a sequence to be predicted
            peptides: a list of arbitrary peptides instead of single sequence
           Returns:
            pandas dataframe
           """

        if sequence is not None:
            tempf = os.path.join(self.temppath, name+'.fa')
            seqfile = write_fasta(sequence, id=name, filename=tempf)
        elif peptides is not None:
            if type(peptides) is not pd.DataFrame:
                peptides = self.seqs_to_dataframe(peptides)
            #print (peptides[:3])
            tempf = tempfile.mktemp()+'.txt'
            seqfile = sequtils.dataframe_to_fasta(peptides, outfile=tempf, seqkey='peptide', idkey='name')
            length = peptides.peptide.str.len().max()
        else:
            return

        if not os.path.exists(self.iedbmhc1_path):
            print ('IEDB MHC-I tools not found, set the iedbmhc1path variable')
            return
        if method == None:
            method = self.method
        if method not in self.methods:
            print ('available methods are %s' %self.methods)
            return
        self.method = method

        cmd = os.path.join(self.iedbmhc1_path,'src/predict_binding.py')
        cmd = cmd+' %s %s %s %s' %(method,allele,length,seqfile)
        if show_cmd == True:
            print (cmd)
        from subprocess import Popen, PIPE
        try:
            p = Popen(cmd, stdout=PIPE, shell=True)
            temp,error = p.communicate()
            #print (temp)
        except OSError as e:
            print (e)
            return
        #print (temp)
        df = self.prepare_data(temp, name)
        return df

    def prepare_data(self, rows, name):
        """Prepare data from results"""

        try:
            df = pd.read_csv(io.BytesIO(rows),sep="\t")
        except:
            #print ('no output to read')
            return
        #print (df)
        if len(df)==0:
            print (rows) #should print error string from output
            return
        df = df.replace('-',np.nan)
        df = df.dropna(axis=1,how='all')
        df.rename(columns={'start':'pos'},
                           inplace=True)
        df=df.drop(columns=['seq_num','end'])
        df['core'] = df.peptide
        df['name'] = name
        if 'method' not in df.columns:
            df['method'] = self.method
        if self.method in ['IEDB_recommended','consensus']:
            df['ic50'] = df.filter(regex="ic50").mean(1)
        if not 'score' in df.columns:
            df['score'] = df.ic50#.apply( lambda x: 1-math.log(x, 50000))
        self.get_ranking(df)
        self.data = df
        #print (df[:10])
        return df

    def get_alleles(self):
        """Get available alleles from model_list file and
            convert to standard names"""

        alleles = list(pd.read_csv(os.path.join(datadir, 'iedb_mhc1_alleles.csv')).allele.values)
        #alleles = sorted(list(set([get_standard_mhc1(i) for i in alleles])))
        return alleles

    def get_allele_data(self):
        path = self.iedbmhc1_path
        if not os.path.exists(path):
            print ('iedb tools not found')
            return
        c = os.path.join(path,'src/predict_binding.py')
        for m in self.methods:
            cmd = c + ' %s mhc' %m
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
            print (temp)
        return

    def check_install(self):
        path = self.iedbmhc1_path
        if not os.path.exists(path):
            print ('iedb mhcII tools not found')
            return
        cmd = os.path.join(path,'src/predict_binding.py')
        try:
            temp = subprocess.check_output(cmd, shell=True)
            print('IEDB MHC-I tools appears to be installed')
            return 1
        except:
            return 0

    @classmethod
    def _get_name(self):
        return 'iedbmhc1'

class IEDBMHCIIPredictor(Predictor):
    """Using IEDB MHC-II method, requires tools to be installed locally"""

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'iedbmhc2'
        #score will store ic50 value for most methods except sturniolo
        self.scorekey = 'score'
        self.cutoff = 500
        self.operator = '<'
        self.rankascending = 1
        self.methods = ['comblib','consensus3','IEDB_recommended',
                        'NetMHCIIpan','nn_align','smm_align','sturniolo']
        self.method = 'IEDB_recommended'
        self.iedbmhc2_path = defaults['iedbmhc2_path']
        self.qf = self.get_quantile_data()
        return

    def prepare_data(self, rows, name):
        """Read data from raw output"""

        if len(rows) == 0:
            return

        df = pd.read_csv(io.BytesIO(rows),sep=r'\t',engine='python',index_col=False)
        #print (df.columns)
        extracols = ['Start','End','comblib_percentile','smm_percentile','nn_percentile',
                     'Sturniolo core',' Sturniolo score',' Sturniolo percentile']
        df.reset_index(inplace=True)
        df.rename(columns={'index':'pos','Sequence': 'peptide','Allele':'allele'},
                           inplace=True)
        df['core'] = df[df.filter(regex="core").columns[0]]
        df['name'] = name

        if 'method' not in df.columns:
            df['method'] = self.method
        #print (df.loc[0])
        method = df.loc[0].method
        if method == 'Sturniolo':
            df['score'] = df.sturniolo_score
        else:
            df['ic50'] = df.filter(regex="ic50").mean(1)
        if not 'score' in df.columns:
            df['score'] = df.ic50#.apply( lambda x: 1-math.log(x, 50000))

        self.get_ranking(df)
        self.data = df
        return df

    def predict(self, sequence=None, peptides=None, length=15, overlap=None, show_cmd=False,
                   allele='HLA-DRB1*01:01', method='IEDB_recommended', name='', **kwargs):
        """Use IEDB MHC-II python module to get predictions. Requires that the IEDB MHC-II
        tools are installed locally. A sequence argument is provided since the cmd line
        only accepts whole sequence to be fragmented.
        """

        self.method = method
        if sequence is not None:
            seqfile = tempfile.mktemp()+'.fa'
            write_fasta(sequence, id=name, filename=seqfile)
        elif peptides is not None:
            if type(peptides) is not pd.DataFrame:
                peptides = self.seqs_to_dataframe(peptides)
            #print (peptides[:3])
            tempf = tempfile.mktemp()+'.txt'
            seqfile = sequtils.dataframe_to_fasta(peptides, outfile=tempf, seqkey='peptide', idkey='name')
            length = peptides.peptide.str.len().max()
        else:
            return

        path = self.iedbmhc2_path
        if method == None: method = 'IEDB_recommended'
        if not os.path.exists(path):
            print ('iedb MHC-II tools not found')
            return
        cmd = os.path.join(path,'mhc_II_binding.py')
        cmd = cmd+' %s %s %s' %(method,allele,seqfile)
        if show_cmd is True:
            print (cmd)
        try:
            #temp = subprocess.check_output(cmd, shell=True)
            from subprocess import Popen, PIPE, STDOUT
            p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
            temp = p.stdout.read()
            #print (temp)
        except Exception as e:
            print (e)
            print ('check allele %s is correct.' %allele)
            return

        data = self.prepare_data(temp, name)
        return data

    def get_alleles(self):
        alleles = list(pd.read_csv(os.path.join(datadir, 'iedb_mhc2_alleles.csv')).allele.values)
        return alleles

    def check_install(self):
        path = self.iedbmhc2_path
        if not os.path.exists(path):
            print ('iedb mhcII tools not found')
            return
        cmd = os.path.join(path,'mhc_II_binding.py')
        try:
            temp = subprocess.check_output(cmd, shell=True)
            print('IEDB MHC-II tools appears to be installed')
            return 1
        except:
            return 0

    @classmethod
    def _get_name(self):
        return 'iedbmhc2'

class TEpitopePredictor(Predictor):
    """Predictor using TepitopePan QM method"""

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'tepitope'
        self.pssms = tepitope.get_pssms()
        self.cutoff = 2
        self.operator = '>'
        self.rankascending = 0
        #load precalculated quantiles for sample peptides
        self.qf = self.get_quantile_data()

    def predict(self, peptides=None, length=9, overlap=1,
                 allele='HLA-DRB1*0101', name='', pseudosequence=None,
                 **kwargs):

        if length not in [9,11,12,15]:
            print ('you should use lengths of 9 or 11, 13 or 15.')
        allele = allele.replace(':','')
        if not allele in self.pssms:
            #print 'computing virtual matrix for %s' %allele
            m = tepitope.create_virtual_pssm(allele)
            if m is None:
                print ('no such allele', allele)
                return pd.DataFrame()
        else:
            m = self.pssms[allele]
        m = m.transpose().to_dict()
        result = tepitope.get_scores(m, peptides=peptides)
        df = self.prepare_data(result, name, allele)
        self.data = df
        #print(df[:5])
        return df

    def supported_lengths(self):
        """Return supported peptide lengths"""
        return [9,11,13,15]

    def get_alleles(self):
        return tepitope.get_alleles()

    def check_alleles(self, alleles):
        alleles = [i.replace(':','') for i in alleles]
        a = self.get_alleles()
        found = list((set(a) & set(alleles)))
        return found

    @classmethod
    def _get_name(self):
        return 'tepitope'

'''class IEDBBCellPredictor(Predictor):
    """Using IEDB tools methods, requires iedb bcell tools.
       see http://tools.immuneepitope.org/bcell """

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'iedbbcell'
        self.scorekey = 'Score'
        self.methods = ['Chou-Fasman', 'Emini', 'Karplus-Schulz',
                        'Kolaskar-Tongaonkar', 'Parker', 'Bepipred']
        self.cutoff = 0.9
        self.operator = '>'
        self.rankascending = 0
        self.method = 'Bepipred'
        self.path = iedbbcellpath

    def predict(self, sequence=None, peptides=None, window=None, name=''):
        """Uses code from iedb predict_binding.py """

        value = self.method
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
        self.prepare_data(temp, name)
        os.chdir(currpath)
        return self.data

    def prepare_data(self, temp, name):

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

    def predict_sequences(self, recs, names=None, save=False,
                        label='', path='', **kwargs):
        """Get predictions for a set of proteins - no alleles so we override
        the base method for this too. """

        recs = sequtils.get_cds(recs)
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
        return'''

class MHCFlurryPredictor(Predictor):
    """
    Predictor using MHCFlurry for MHC-I predictions. Requires you to
    install the python package mhcflurry with dependencies.
    see https://github.com/hammerlab/mhcflurry
    """

    def __init__(self, data=None, **kwargs):
        Predictor.__init__(self, data=data)
        self.name = 'mhcflurry'
        self.cutoff = 500
        self.operator = '<'
        self.scorekey = 'score'
        self.rankascending = 1
        if self.check_install() == False:
            return
        self._check_models()
        from mhcflurry import Class1AffinityPredictor
        self.predictor = Class1AffinityPredictor.load()
        #load precalculated quantiles for sample peptides
        self.qf = self.get_quantile_data()
        return

    def predict(self, peptides=None, overlap=1,
                      allele='HLA-A0101', name='', **kwargs):
        """Uses mhcflurry python classes for prediction"""

        try:
            df = self.predictor.predict_to_dataframe(peptides=peptides, allele=allele)
        except:
            print ('failed to predict for allele %s' %allele)
            return
        #print (df[:5])
        df = self.prepare_data(df, name, allele)
        self.data = df
        return df

    def prepare_data(self, df, name, allele):
        """Post process dataframe to alter some column names"""

        df = df.rename(columns={'Allele':'allele','Peptide':'peptide'})
        df['name'] = name
        df['pos'] = df.index
        df['score'] = df.prediction
        df = df.drop(columns=['prediction_low', 'prediction_high'])
        self.get_ranking(df)
        return df

    def convert_allele_name(self, r):
        if ':' not in r:
            return r[:5]+'*'+r[5:7]+':'+r[7:]
        else:
            return r[:5]+'*'+r[5:7]+r[7:]

    def get_alleles(self):
        with open(os.path.join(datadir,'mhcflurry_alleles.txt')) as f:
            p = f.readlines()
        p = [x.strip() for x in p]
        p = list(filter(None, p))
        return p

    def check_install(self):
        try:
            from mhcflurry import Class1AffinityPredictor
            return True
        except:
            #print ('use pip install mhcflurry to install.')
            return False

    def _check_models(self):
        try:
            import mhcflurry
            mhcflurry.class1_affinity_predictor.get_default_class1_models_dir()
        except:
            print ('first use. downloading class-I models')
            subprocess.check_output('mhcflurry-downloads fetch models_class1', shell=True)
            return True
        return False

    @classmethod
    def _get_name(self):
        return 'mhcflurry'

class MHCNuggetsPredictor(Predictor):
    """
    Predictor using MHCNuggets for MHC-I predictions. Requires you to
    install the package locally from https://github.com/KarchinLab/mhcnuggets
    see https://www.biorxiv.org/content/early/2017/07/27/154757
    """

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'mhcnuggets'
        self.cutoff = 500
        self.operator = '<'
        self.scorekey = 'ic50'
        self.rankascending = 1
        self.qf = self.get_quantile_data()
        return

    def predict(self, peptides=None, overlap=1,
                      allele='HLA-A01:01', name='temp', **kwargs):
        """Uses cmd line call to mhcnuggets."""

        allele = allele.replace('*','')
        from mhcnuggets.src.predict import predict
        pepfile = self.write_seqs(peptides)
        outfile = tempfile.mktemp()+'.txt'
        #print(outfile)
        predict(class_='I', peptides_path=pepfile, mhc=allele, output=outfile)
        res = pd.read_csv(outfile)
        df = self.prepare_data(res, name, allele)
        return df

    def prepare_data(self, df, name, allele):
        """Get result into dataframe"""

        df['name'] = name
        df['allele'] = allele
        df['pos'] = df.index.copy()
        self.get_ranking(df)
        self.data = df
        return df

    def write_seqs(self, peptides):
        tempf = tempfile.mktemp()+'.pep'
        f = open(tempf,'w')
        for p in peptides:
            f.write("%s\n" % p)
        f.close()
        return tempf

    def get_alleles(self):
        with open(os.path.join(datadir,'mhcnuggets_alleles.txt')) as f:
            p = f.readlines()
        p = [x.strip() for x in p]
        p = list(filter(None, p))
        return p

    @classmethod
    def _get_name(self):
        return 'mhcnuggets'

class DummyPredictor(Predictor):
    """Returns random scores. Used for testing"""
    def __init__(self, data=None, scoring=None):
        Predictor.__init__(self, data=data)
        self.name = 'null'
        self.scorekey = 'score'
        self.cutoff = 500
        self.operator = '<'
        self.rankascending = 1

    def predict(self, peptides, allele='HLA-A*01:01', name='temp',
             **kwargs):

        sc = np.random.randint(0,1,len(peptides))
        res = pd.DataFrame(np.column_stack([peptides,sc]),columns=['peptide','log50k'])
        res['score'] = res.log50k.astype('float').apply(lambda x: log50k2aff(x))
        df = self.prepare_data(res,name,allele)
        return df

    @classmethod
    def _get_name(self):
        return 'dummy'

class BasicMHCIPredictor(Predictor):
    """Built-in basic MHC-I predictor. Should be used as a fallback if no other
       predictors available."""

    def __init__(self, data=None, scoring=None):
        Predictor.__init__(self, data=data)
        self.name = 'basicmhc1'
        self.scorekey = 'score'
        self.cutoff = 500
        self.operator = '<'
        self.rankascending = 1
        self.qf = self.get_quantile_data()
        return

    def predict(self, peptides, allele='HLA-A*01:01', name='temp', **kwargs):
        """Encode and predict peptides with daved regressor"""

        encoder = mhclearn.blosum_encode
        if type(peptides) is list:
            s = pd.Series(peptides)
        else:
            s = peptides
        #encode
        X = s.apply(lambda x: pd.Series(encoder(x)),1)
        reg = mhclearn.get_model(allele)
        if reg is None:
            print ('no model for this allele')
            return
        sc = reg.predict(X).round(4)
        res = pd.DataFrame(np.column_stack([peptides,sc]),columns=['peptide','log50k'])
        res['score'] = res.log50k.astype('float').apply(lambda x: mhclearn.log50k2aff(x)).round(2)
        #print(res.dtypes)
        df = self.prepare_data(res,name,allele)
        return df

    def prepare_data(self, df, name, allele):
        """Put raw prediction data into DataFrame and rank,
           override for custom processing. Can be overriden for
           custom data."""

        df['name'] = name
        df['allele'] = allele
        self.get_ranking(df)
        return df

    def supported_lengths(self):
        """Return supported peptide lengths"""
        return [9]

    def get_alleles(self):
        p = mhclearn.get_allele_names()
        return p

    def check_install(self):

        reg = mhclearn.get_predictor('HLA-A*01:01')
        if reg is None:
            return False
        return True

    @classmethod
    def _get_name(self):
        return 'basicmhc1'
