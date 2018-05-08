#!/usr/bin/env python

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

"""
    epitopepredict config
    Created March 2016
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, string, time
import types, re, subprocess, glob, shutil
from collections import OrderedDict
import pandas as pd
try:
    import configparser
except:
    import ConfigParser as configparser

path = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(path, 'data')
home = os.path.expanduser("~")
config_path = os.path.join(home, '.epitopepredict')

baseoptions = OrderedDict()
baseoptions['base'] = {'predictors': 'tepitope',
                'mhc2_alleles':'HLA-DRB1*01:01,HLA-DRB1*04:01',
                'mhc1_alleles':'HLA-A*01:01',
                'mhc1_length': 11,
                'mhc2_length': 11,
                'n': 2, #number of alleles
                'cutoff_method': 'default',
                'cutoffs': .95, #percentile cutoff
                'sequence_file':'', #genbank/fasta file
                'peptide_file':'', #plain text list of peptides
                'path': 'results',
                'overwrite': 'no',
                'verbose':'no',
                'names': '', #subset of protein names
                'overwrite': 'no',
                'cpus': 1,
                'compression': '',
                'fasta_header_sep': ' '}

baseoptions['iedbtools'] = {'iedbmhc1_path':'', 'iedbmhc2_path':'',
                            'iedb_mhc1_method':'IEDB_recommended',
                            'iedb_mhc2_method':'IEDB_recommended'}

baseoptions['neopredict'] = {'vcf_files':'','maf_files':''}

def write_default_config():
    """Write a default config to users .config folder. Used to add global settings."""

    fname = os.path.join(config_path, 'default.conf')
    if not os.path.exists(fname):
        try:
            #os.mkdir(config_path)
            os.makedirs(config_path)
        except:
            pass
        write_config(conffile=fname, defaults=baseoptions)
    return fname

def write_config(conffile='default.conf', defaults={}):
    """Write a default config file"""

    if not os.path.exists(conffile):
        cp = create_config_parser_from_dict(defaults, ['base','iedbtools'])
        cp.write(open(conffile,'w'))
        print ('wrote config file %s' %conffile)
    return conffile

def create_config_parser_from_dict(data=None, sections=['base','iedbtools'], **kwargs):
    """Helper method to create a ConfigParser from a dict of the form shown in
       baseoptions"""

    if data is None:
        data = baseoptions
    #print (data)
    cp = configparser.ConfigParser()
    for s in sections:
        cp.add_section(s)
        if not data.has_key(s):
            continue
        for name in sorted(data[s]):
            val = data[s][name]
            if type(val) is list:
                val = ','.join(val)
            cp.set(s, name, str(val))

    #use kwargs to create specific settings in the appropriate section
    for s in cp.sections():
        opts = cp.options(s)
        for k in kwargs:
            if k in opts:
                cp.set(s, k, kwargs[k])
    return cp

def parse_config(conffile=None):
    """Parse a configparser file"""

    f = open(conffile,'r')
    cp = configparser.ConfigParser()
    try:
        cp.read(conffile)
    except Exception as e:
        print ('failed to read config file! check format')
        print ('Error returned:', e)
        return
    f.close()
    return cp

def get_options(cp):
    """Makes sure boolean opts are parsed"""

    from collections import OrderedDict
    options = OrderedDict()
    #options = cp._sections['base']
    for section in cp.sections():
        options.update( (cp._sections[section]) )
    for o in options:
        for section in cp.sections():
            try:
                options[o] = cp.getboolean(section, o)
            except:
                pass
            try:
                options[o] = cp.getint(section, o)
            except:
                pass
    return options

def print_options(options):
    """Print option key/value pairs"""

    for key in options:
        print (key, ':', options[key])
    print ()

def check_options(opts):
    """Check for missing default options in dict. Meant to handle
       incomplete config files"""

    sections = list(baseoptions.keys())
    for s in sections:
        defaults = dict(baseoptions[s])
        for i in defaults:
            if i not in opts:
                opts[i] = defaults[i]
    return opts
