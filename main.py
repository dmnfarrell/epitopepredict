#!/usr/bin/env python

"""
    MHC prediction command line script
    Created March 2016
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os
import pandas as pd
from configparser import ConfigParser
from collections import OrderedDict
import epitopepredict as ep
from epitopepredict import base, analysis, sequtils, tests

defaultpath = os.getcwd()
optvalues = (('predictors', 'tepitope'),
               ('mhc2alleles','HLA-DRB1*01:01,HLA-DRB1*04:01'),
               ('mhc1alleles','HLA-A*01:01'),
               ('n', 2),
               ('cutoff_method', 'default'),
               ('cutoff',0.98),
               ('genome', 'name.gb'),
               ('path', os.getcwd()),
               ('prefix', '_'), #prefix for subfolders
               ('overwrite', 'no'),
               ('names', ''),
               ('genome_analysis', 'no'))
defaultopts = OrderedDict(optvalues)

def createConfigParserfromOptions(opts, section):
    """Helper method to create a ConfigParser from a dict of options"""

    cp = ConfigParser()
    s='prediction'
    cp.add_section(s)
    print('writing a new config file')
    for name in opts:
        val = opts[name]
        print(name,val)
        cp.set(s, name, str(val))
    #cp.write(open(filename,'w'))
    return cp

def createConfig(opts=None, conffile='default.conf'):
    """Create a basic config file with default options and/or custom values"""

    if opts == None:
        opts = defaultopts
    c = ConfigParser()
    wdir = os.path.join(defaultpath, 'workingdir')
    cp = createConfigParserfromOptions(opts, 'default')
    cp.write(open(conffile,'w'))

    #self.parseConfig(conffile)
    return cp

def parseConfig(conffile=None):
    """Parse the config file"""

    f = open(conffile,'r')
    cp = ConfigParser()
    try:
        cp.read(conffile)
    except:
        pass
    print (cp)
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
    print (data)
    return data

def run(predictors=[], cutoff=0.98, cutoff_method='default',
         mhc2alleles=[], mhc1alleles=[],
         n=2,  genome='',
         path='', prefix='_',
         overwrite=False,
         genome_analysis=False,
         names = ''):
    """Run a workflow using config settings"""

    genome = sequtils.genbank2Dataframe(genome, cds=True)
    predictors = predictors.split(',')
    mhc2alleles = mhc2alleles.split(',')
    print (mhc2alleles)
    for p in predictors:
        print (p)
        P = ep.getPredictor(p)
        savepath = os.path.join(path,prefix+p)
        P.predictProteins(genome, length=11, alleles=mhc2alleles, #names=names,
                          path=savepath, overwrite=overwrite)
        P.load(path=savepath)
        #pb = P.getPromiscuousBinders(n=n, perc=cutoff, cutoff_method=cutoff_method)

    return

def main():
    "Run the application"
    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-c", "--config", dest="config",
                        help="Configuration file", metavar="FILE")
    parser.add_option("-r", "--run", dest="run",  action="store_true",
                        default=False, help="Run")
    #parser.add_option("-t", "--test", dest="test",  action="store_true",
    #                    default=False, help="tests")
    opts, remainder = parser.parse_args()
    if opts.config != None:
        #from epitopepredict.tests import *
        #unittest.main()
        cp = parseConfig(opts.config)
        #print (cp)
    else:
        createConfig()
    if opts.run == True:
        kwargs = config2Dict(cp)['prediction']
        run(**kwargs)


if __name__ == '__main__':
    main()
