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
blosum=pd.read_csv(os.path.join(datadir,'blosum62.csv'))

def aff2log50k(a):
    return 1 - (math.log(a) / math.log(50000))

def log50k2aff(a):
    return np.power(50000,1-a)

def blosum_encode(seq):
    s=list(seq)
    x = pd.DataFrame([blosum[i] for i in seq]).reset_index(drop=True)
    #show_matrix(x)
    e = x.values.flatten()
    return e

def nlf_encode(seq):
    x = pd.DataFrame([nlf[i] for i in seq]).reset_index(drop=True)
    #print (x)
    e = x.values.flatten()
    return e

def simple_score(true,sc,cutoff=None):

    if cutoff!=None:
        true = (true<=cutoff).astype(int)
        sc = (sc<=cutoff).astype(int)
        #print (sc.value_counts())
    fpr, tpr, thresholds = metrics.roc_curve(true, sc, pos_label=1)
    r = metrics.auc(fpr, tpr)
    #print (r)
    return  r

def get_protein_set():

    syf = os.path.join(datadir, 'SYF_set.fasta')
    return sequtils.fasta_to_dataframe(syf)

def get_training_set(allele=None):

    b = pd.read_csv(os.path.join(datadir, 'curated_training_data.no_mass_spec.zip'))
    df=b
    if allele is not None:
        df = b.loc[b.allele==allele].copy()
    #print (len(b),len(df))
    df['log50k'] = df.ic50.apply(lambda x: aff2log50k(x))
    df['length'] = df.peptide.str.len()
    df = df[(df.length==9) & (df.ic50<50000) ]
    df['binder'] = (df.ic50<500).astype(int)
    return df

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

def get_predictor(allele):
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
