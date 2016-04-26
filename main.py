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
from epitopepredict import base, analysis, sequtils, tests

defaultpath = os.getcwd()
optvalues = (('predictors', 'tepitope'),
               ('mhc2alleles','HLA-DRB1*01:01,HLA-DRB1*04:01'),
               ('mhc1alleles','HLA-A*01:01'),
               ('n', 2), ('cutoff_method', 'default'),
               ('genome', 'name.gb'),
               ('path', os.getcwd()),
               ('find_clusters',True),
               ('genome_analysis', False))
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
    return

def run():
    """Run a workflow using config settings"""

    genome = sequtils.genbank2Dataframe(genomefile, cds=True)
    P1 = ep.getPredictor('tepitope')
    return


def main():
    "Run the application"
    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-c", "--config", dest="config",
                        help="Configuration file", metavar="FILE")
    #parser.add_option("-t", "--test", dest="test",  action="store_true",
    #                    default=False, help="Run tests")
    opts, remainder = parser.parse_args()
    if opts.config != None:
        #from epitopepredict.tests import *
        #unittest.main()
        cp = parseConfig(opts.config)
        print (cp)
    else:
        createConfig()


if __name__ == '__main__':
    main()
