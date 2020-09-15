#!/usr/bin/env python

"""
    Basic MHC binding prediction predictors
    Created October 2018
    Copyright (C) Damien Farrell
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

from __future__ import absolute_import, print_function
import sys, os, math
import numpy as np
import pandas as pd
from collections import OrderedDict
from epitopepredict import peptutils, sequtils

home = os.path.expanduser("~")
config_path = os.path.join(home, '.config/epitopepredict')
module_path = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(module_path, 'mhcdata')
models_path = os.path.join(config_path, 'models')
nlf = pd.read_csv(os.path.join(datadir,'NLF.csv'),index_col=0)
blosum62 = pd.read_csv(os.path.join(datadir,'blosum62.csv'))
blosum50 = pd.read_csv(os.path.join(datadir,'blosum50.csv'))
codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
         'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def check_valid_prot_seq(seq):
    for i in seq:
        if i not in mhclearn.codes:
            return False

def aff2log50k(a):
    return 1 - (math.log(a) / math.log(50000))

def log50k2aff(a):
    return np.power(50000,1-a)

def random_encode(p):
    """Random encoding of a peptide for testing"""

    return [np.random.randint(20) for i in pep]

def one_hot_encode(seq):
    """One hot encoding of peptide"""

    o = list(set(codes) - set(seq))
    s = pd.DataFrame(list(seq))
    x = pd.DataFrame(np.zeros((len(seq),len(o)),dtype=int),columns=o)
    a = s[0].str.get_dummies(sep=',')
    a = a.join(x)
    a = a.sort_index(axis=1)
    #a = a.set_index([a.index,list(seq)])
    #show_matrix(a)
    e = a.values.flatten()
    return e

def blosum_encode(seq):
    """Blosum62 encoding of peptide"""

    s=list(seq)
    x = pd.DataFrame([blosum62[i] for i in seq]).reset_index(drop=True)
    x = x.iloc[:,:-4]
    e = x.values.flatten()
    return e

def nlf_encode(seq):
    """NLF encoding of a peptide (from Nanni and Lumini)"""

    x = pd.DataFrame([nlf[i] for i in seq]).reset_index(drop=True)
    e = x.values.flatten()
    return e

def auc_score(true, sc, cutoff=None):

    from sklearn import metrics
    if cutoff!=None:
        true = (true<=cutoff).astype(int)
        sc = (sc<=cutoff).astype(int)
    fpr, tpr, thresholds = metrics.roc_curve(true, sc, pos_label=1)
    r = metrics.auc(fpr, tpr)
    return  r

def get_protein_set():

    syf = os.path.join(datadir, 'SYF_set.fasta')
    return sequtils.fasta_to_dataframe(syf)

def get_training_set(allele=None, length=None):
    """Get training set for MHC-I data."""

    b = pd.read_csv(os.path.join(datadir, 'curated_training_data.no_mass_spec.zip'))
    eval1 = get_evaluation_set1()
    df = b.loc[~b.peptide.isin(eval1.peptide)].copy()
    if allele is not None:
        df = b.loc[b.allele==allele].copy()

    df['log50k'] = df.ic50.apply(lambda x: aff2log50k(x))
    df['length'] = df.peptide.str.len()
    if length != None:
        df = df[(df.length==length)]
    df = df[df.ic50<50000]
    df = df[df.measurement_type=='quantitative']
    #df['binder'] = df.loc[df.ic50<500].astype(int)
    return df

def get_evaluation_set1(allele=None, length=None):
    """Get eval set of peptides"""

    e = pd.read_csv(os.path.join(datadir, 'binding_data_2013.zip'),comment='#')
    if allele is not None:
        e = e[e.allele==allele]
    if length != None:
        e = e[(e.length==length) ]
    e['log50k'] = e.ic50.apply(lambda x: aff2log50k(x)).round(2)
    return e

def get_allele_names():
    """Get allele names from training set"""

    b = get_training_set()
    a = b.allele.value_counts()
    a = a[a>100]
    return list(a.index)

def build_predictor(allele, length=9, encoder=None):
    """Build a regression model using peptide encoder and test data"""

    from sklearn.neural_network import MLPRegressor
    #get allele specific predictor
    if encoder==None:
        encoder=one_hot_encode
    data = get_training_set(allele, length=length)
    if len(data)<100:
        return
    reg = MLPRegressor(hidden_layer_sizes=(1,), alpha=0.0005, max_iter=1500,
                        activation='logistic', solver='lbfgs', random_state=2)
    #reg = LinearRegression()
    X = data.peptide.apply(lambda x: pd.Series(encoder(x)),1)
    y = data.log50k
    print ('trained model:', allele, len(X), length)
    reg.fit(X,y)
    return reg

def train_models(overwrite=False, alleles=None, encoder=None):
    """Train and save models. May be needed if version of joblib is different."""

    os.makedirs(models_path, exist_ok=True)
    #encoders = {'blosum':blosum_encode, 'onehot':one_hot_encode, 'nlf':nlf_encode}
    import joblib
    if alleles == None:
        alleles = get_allele_names()
    if type(alleles) is str:
        alleles = [alleles]
    if encoder == None:
        encoder = one_hot_encode
    #else:
    #    encoder = encoders[encoder]

    lengths = [9,10,11]
    for a in alleles:
        for l in lengths:
            fname = os.path.join(models_path, a+'_'+str(l)+'.joblib')
            #print (fname)
            if os.path.exists(fname) and overwrite == False:
                continue
            reg = build_predictor(a, l, encoder)
            if reg is not None:
                joblib.dump(reg, fname, protocol=2)
    return

def get_model(allele, length):
    """Get a regression model."""

    try:
        import sklearn
    except:
        print ('you need scikit-learn to use this predictor')
    import joblib
    #os.makedirs(models_path, exist_ok=True)
    fname = os.path.join(models_path, allele+'_'+str(length)+'.joblib')
    if os.path.exists(fname):
        reg = joblib.load(fname)
        return reg
