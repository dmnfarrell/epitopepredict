#!/usr/bin/env python

"""
    Basic mhc binding prediction predictors
    Created October 2018
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, math
import numpy as np
import pandas as pd
from collections import OrderedDict
from epitopepredict import peptutils, sequtils

module_path = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(module_path, 'mhcdata')
modelspath = os.path.join(module_path, 'models')
nlf = pd.read_csv(os.path.join(datadir,'NLF.csv'),index_col=0)
blosum62 = pd.read_csv(os.path.join(datadir,'blosum62.csv'))

def aff2log50k(a):
    return 1 - (math.log(a) / math.log(50000))

def log50k2aff(a):
    return np.power(50000,1-a)

def random_encode(p):
    return [np.random.randint(20) for i in pep]

def blosum_encode(seq):
    s=list(seq)
    x = pd.DataFrame([blosum62[i] for i in seq]).reset_index(drop=True)
    #show_matrix(x)
    e = x.values.flatten()
    return e

def nlf_encode(seq):
    x = pd.DataFrame([nlf[i] for i in seq]).reset_index(drop=True)
    #print (x)
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
    eval1 = get_evaluation_set()
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

def get_evaluation_set(allele=None, length=None):
    """Get eval set of peptides"""

    e = pd.read_csv(os.path.join(datadir, 'eval_set1.csv'))
    #remove evaluation peptides
    if allele is not None:
        e = e[e.allele==allele]
    e['length'] = e.peptide.str.len()
    if length != None:
        e = e[(e.length==length) ]
    e['ic50'] = e.log50k.apply(log50k2aff)
    return e

def get_allele_names():

    b = get_training_set()
    a = b.allele.value_counts()
    a = a[a>200]
    return list(a.index)

def build_predictor(allele):
    """Build a regression model using peptide encoder and test data"""

    #get allele specific predictor
    encoder=nlf_encode
    #encoder=pc19_encode
    data = get_test_data(allele)
    from sklearn.linear_model import LinearRegression
    reg = LinearRegression()
    X = data.sequence.apply(lambda x: pd.Series(encoder(x)),1)
    y = data.log50k
    print (allele, len(X))
    reg.fit(X,y)
    return reg

def get_model(allele):
    """Get a regression model."""

    try:
        import sklearn
    except:
        print ('you need scikit-learn to use this predictor')
    from sklearn.externals import joblib
    fname = os.path.join(modelspath, allele+'.joblib')
    if os.path.exists(fname):
        reg = joblib.load(fname)
        return reg
    else:
        return
