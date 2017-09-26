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
from collections import OrderedDict
from . import base, config, analysis, sequtils, plotting, tests

defaultpath = os.getcwd()

class WorkFlow(object):
    """Class for implementing a rna/mirna workflow from a set of options"""
    def __init__(self, opts={}):
        for i in opts:
            self.__dict__[i] = opts[i]
        return

    def setup(self):
        """Setup main parameters"""

        base.iedbmhc1path = self.iedbmhc1_path
        base.iedbmhc2path = self.iedbmhc2_path
        self.sequences = get_sequences(self.sequence_file)
        self.mhc1_alleles = self.mhc1_alleles.split(',')
        self.mhc2_alleles = self.mhc2_alleles.split(',')
        if len(self.mhc1_alleles)==0 and len(self.mhc2_alleles)==0:
            return False
        self.predictors = self.predictors.split(',')
        for p in self.predictors:
            if p not in base.predictors:
                print ('unknown predictor in config file. Use:')
                show_predictors()
                return False
        if self.mhc1_alleles[0] in base.mhc1_presets:
            self.mhc1_alleles = base.get_preset_alleles(self.mhc1_alleles[0])
        elif self.mhc2_alleles[0] in base.mhc2_presets:
            self.mhc2_alleles = base.get_preset_alleles(self.mhc2_alleles[0])

        self.cutoff = float(self.cutoff)
        self.names = self.names.split(',')
        if self.names == ['']: self.names=None
        if not os.path.exists(self.path) and self.path != '':
            os.mkdir(self.path)
        return True

    def run(self):
        """Run workflow"""

        preds = []
        for p in self.predictors:
            P = base.get_predictor(p)
            preds.append(P)
            savepath = os.path.join(self.path, p)
            if self.overwrite == True and os.path.exists(savepath):
                shutil.rmtree(savepath)
            if p in ['iedbmhc1','mhcflurry']:
                a = self.mhc1_alleles
                length = self.mhc1_length
                check_mhc1_length(length)
                method = self.iedb_mhc1_method
            else:
                a = self.mhc2_alleles
                length = self.mhc2_length
                method = self.iedb_mhc2_method
            if method == '': method = None
            print ('predictor:', p, method)
            print ('alleles:',a)
            if p == 'iedbmhc1' and check_iedbmhc1_path() == False:
                continue

            P.predictProteins(self.sequences, length=length, alleles=a, names=self.names,
                              path=savepath, overwrite=self.overwrite, verbose=self.verbose,
                              method=method)
            #load results into predictor
            P.load(path=savepath)
            if P.data is None:
                print ('no results were found, did predictor run?')
                return
            cutoff = self.cutoff
            cutoff_method = self.cutoff_method
            n = self.n
            b = P.getBinders(cutoff=cutoff, value=cutoff_method)
            b.to_csv(os.path.join(self.path,'binders_%s_%s.csv' %(p,n)))

            pb = P.promiscuousBinders(n=int(n), cutoff=cutoff, value=cutoff_method)
            print ('found %s promiscuous binders at cutoff %s' %(len(pb),cutoff))
            pb.to_csv(os.path.join(self.path,'prom_binders_%s_%s.csv' %(p,n)))
            if self.verbose == True:
                print ('top promiscuous binders:')
                print (pb[:10])
            if self.genome_analysis == True:
                cl = analysis.find_clusters(pb, genome=self.sequences)
                cl.to_csv(os.path.join(self.path,'clusters_%s.csv' %p))
            print ('-----------------------------')

        self.preds = preds
        if self.plots == True:
            self.plot_results()
        return

    def plot_results(self):
        """Plot results of predictions"""

        preds = self.preds
        prots = self.sequences.locus_tag
        import pylab as plt
        height = 2*len(preds)
        for prot in prots:
            ax = plotting.plot_tracks(preds,name=prot,n=2,cutoff=self.cutoff,
                                          figsize=(14,height),legend=True)
            #plotting.mpl_plot_regions(coords, ax, color='gray')
            #ax = plotting.plot_bars(preds[0], prot, cutoff=20, chunks=1)
            plt.tight_layout()
            plt.savefig('plots/%s.png'%prot, dpi=150)
        print ('saved plots')
        return

    def analysis(self, path):

        prots = self.sequences.locus_tag
        for p in base.predictors:
            P = base.get_predictor(p)
            P.load(os.path.join(path, p))
            print (P)
        return

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

def check_mhc1_length(l):
    if l<9 or l>13:
        print ('use MHCI n-mer lengths from 9-13')
        return False

def check_iedbmhc1_path():
    if not os.path.exists(base.iedbmhc1path):
        print ('IEDB MHC tools not found, check path')
        return False

def show_preset_alleles():
    print ('preset allele list ids:')
    for i in base.mhc1_presets+base.mhc2_presets:
        print (i, len( base.get_preset_alleles(i)))

def show_predictors():
    for p in base.predictors:
        print(p)

def print_help():
    print ("""use -h to get options""")

def list_alleles():
    for p in base.predictors:
        print (p)
        print ('-----------------------------')
        P = base.get_predictor(p)
        x = P.getAlleles()
        if type(x) is list:
            for i in x: print (i)
        print ()
    return

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
    parser.add_option("-l", "--list-alleles", dest="list_alleles",  action="store_true",
                        default=False, help="List available alleles")
    parser.add_option("-t", "--test", dest="test",  action="store_true",
                        default=False, help="Do quick test")
    parser.add_option("-a", "--analysis", dest="analysis",
                        help="Analysis path", metavar="FILE")
    opts, remainder = parser.parse_args()

    if opts.config != None:
        cp = config.parse_config(opts.config)
        options = config.get_options(cp)
        options = config.check_options(options)
    else:
        conffile = 'default.conf'
        if not os.path.exists(conffile):
            config.write_default_config(conffile, defaults=config.baseoptions)
    if opts.presets == True:
        show_preset_alleles()
    elif opts.list_alleles == True:
        list_alleles()
    elif opts.run == True:
        W = WorkFlow(options)
        st = W.setup()
        if st == True:
            W.run()
    elif opts.analysis is not None:
        W = WorkFlow()
        W.analysis(opts.analysis)
    else:
        print_help()

if __name__ == '__main__':
    main()
