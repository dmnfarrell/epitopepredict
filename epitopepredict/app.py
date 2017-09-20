#!/usr/bin/env python

"""
    MHC prediction command line script
    Created March 2016
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os
import shutil
import pandas as pd
try:
    import configparser
except:
    import ConfigParser as configparser
from collections import OrderedDict
from . import base, analysis, sequtils, plotting, tests

defaultpath = os.getcwd()
#default configuration values
optvalues = (('predictors', 'tepitope'),
               ('mhc2_alleles','HLA-DRB1*01:01,HLA-DRB1*04:01'),
               ('mhc1_alleles','HLA-A*01:01'),
               ('preset_alleles',''),
               ('mhc1_length', 11),
               ('mhc2_length', 15),
               ('n', 2), #number of alleles
               ('cutoff_method', 'default'),
               ('cutoff',4), #percentile cutoff
               ('sequence_file', ''), #genbank/fasta file
               ('path', 'results'),
               ('overwrite', 'no'),
               ('verbose','no'),
               ('names', ''), #subset of protein names from genome file
               ('overwrite', 'no'),
               ('plots','no'), #whether to save plots
               ('genome_analysis', 'no'))
defaultopts = OrderedDict(optvalues)

def createConfigParserfromOptions(opts, section):
    """Helper method to create a ConfigParser from a dict of options"""

    cp = configparser.ConfigParser()
    s = 'base'
    cp.add_section(s)
    print('writing a new config file')
    for name in opts:
        val = opts[name]
        #print(name,val)
        cp.set(s, name, str(val))
    #cp.write(open(filename,'w'))
    return cp

def createConfig(opts=None, conffile='default.conf'):
    """Create a basic config file with default options and/or custom values"""

    if opts == None:
        opts = defaultopts
    c = configparser.ConfigParser()
    wdir = os.path.join(defaultpath, 'workingdir')
    cp = createConfigParserfromOptions(opts, 'default')
    cp.write(open(conffile,'w'))
    #self.parseConfig(conffile)
    return cp

def parseConfig(conffile=None):
    """Parse the config file"""

    f = open(conffile,'r')
    cp = configparser.ConfigParser()
    try:
        cp.read(conffile)
    except:
        pass
    return cp

def config2Dict(config):
    """Convert confiparser sections to dict"""

    data={}
    for s in config.sections():
        #print (s)
        d = config.items(s)
        data[s]={}
        for i,name in config.items(s):
            #print(i)
            try:
                data[s][i] = (config.getboolean(s,i))
            except:
                data[s][i] = config.get(s,i)
    #print (data)
    return data

def get_sequences(filename):
    """Determine file type and get sequences"""

    ext = os.path.splitext(filename)[1]
    if ext in ['.fa','.faa','.fasta']:
        seqs = sequtils.fasta_to_dataframe(filename)
        print ('found fasta file')
    elif ext in ['.gb','.gbk','.genbank']:
        seqs = sequtils.genbank_to_dataframe(filename, cds=True)
        print ('found genbank file')
    return seqs

def run(predictors=[], cutoff=0.98, cutoff_method='default',
         mhc2_alleles='', mhc1_alleles='', preset_alleles='',
         mhc1_length=11, mhc2_length=15,
         n=2,  sequence_file='',
         path='', #prefix='results',
         overwrite=False,
         plots=False,
         genome_analysis=False,
         names = ''):
    """Run the prediction workflow using config settings"""

    sequences = get_sequences(sequence_file)
    #print (sequences)
    #process these in config2Dict
    predictors = predictors.split(',')
    if preset_alleles in base.mhc1_presets:
        mhc1_alleles = base.get_preset_alleles(preset_alleles)
    elif preset_alleles in base.mhc2_presets:
        mhc2_alleles = base.get_preset_alleles(preset_alleles)
    else:
        mhc1_alleles = mhc1_alleles.split(',')
        mhc2_alleles = mhc2_alleles.split(',')
    cutoff = float(cutoff)

    if not os.path.exists(path):
        os.mkdir(path)
    preds = []
    for p in predictors:
        print ('predictor:', p)
        P = base.get_predictor(p)
        preds.append(P)
        savepath = os.path.join(path, p)
        if overwrite == True and os.path.exists(savepath):
            shutil.rmtree(savepath)
        if p in ['iedbmhc1','mhcflurry']:
            a = mhc1_alleles
            length = mhc1_length
        else:
            a = mhc2_alleles
            length = mhc2_length
        P.predictProteins(sequences, length=length, alleles=a, #names=names,
                          path=savepath, overwrite=overwrite)
        #load into predictor
        P.load(path=savepath)
        b = P.getBinders(cutoff=cutoff)#, value=cutoff_method)
        b.to_csv(os.path.join(path,'binders_%s_%s.csv' %(p,n)))

        pb = P.promiscuousBinders(n=int(n), cutoff=cutoff)
        print ('found %s promiscuous binders at cutoff %s' %(len(pb),cutoff))
        pb.to_csv(os.path.join(path,'prom_binders_%s_%s.csv' %(p,n)))
        if genome_analysis == True:
            cl = analysis.find_clusters(pb, genome=sequences)
            cl.to_csv(os.path.join(path,'clusters_%s.csv' %p))
        print ('-----------------------------')

    #various choices here - we could generate a notebook with the plots
    #embedded ? better than saving all to disk
    prots = sequences.locus_tag
    if plots == True:
        import pylab as plt
        height = 2*len(preds)
        for prot in prots:
            ax = plotting.plot_tracks(preds,name=prot,n=2,cutoff=cutoff,
                                          figsize=(14,height),legend=True)
            #plotting.mpl_plot_regions(coords, ax, color='gray')
            plt.tight_layout()
            plt.savefig('plots/%s.png'%prot, dpi=150)
        print ('saved plots')

    return

def show_preset_alleles():
    print ('preset allele list ids:')
    for i in base.mhc1_presets+base.mhc2_presets:
        print (i, len( base.get_preset_alleles(i)))

def main():
    "Run the application"

    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-c", "--config", dest="config",
                        help="Configuration file", metavar="FILE")
    parser.add_option("-r", "--run", dest="run",  action="store_true",
                        default=False, help="Run the predictions")
    parser.add_option("-p", "--presets", dest="presets",  action="store_true",
                        default=False, help="Show preset allele lists")
    parser.add_option("-t", "--test", dest="test",  action="store_true",
                        default=False, help="Do quick test")
    opts, remainder = parser.parse_args()
    if opts.config != None:
        cp = parseConfig(opts.config)
        #print (cp)
    else:
        if not os.path.exists('default.conf'):
            createConfig()
    if opts.presets == True:
        show_preset_alleles()
    elif opts.run == True:
        kwargs = config2Dict(cp)['base']
        run(**kwargs)

if __name__ == '__main__':
    main()
