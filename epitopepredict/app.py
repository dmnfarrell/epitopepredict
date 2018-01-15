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
from epitopepredict import base, config, analysis, sequtils, plotting, tests

defaultpath = os.getcwd()

class WorkFlow(object):
    """Class for implementing a rna/mirna workflow from a set of options"""
    def __init__(self, opts={}):
        for i in opts:
            self.__dict__[i] = opts[i]
        return

    def setup(self):
        """Setup main parameters"""

        pd.set_option('display.width', 120)
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
        if self.names == None:
            self.names = self.sequences.locus_tag
        for p in self.predictors:
            P = base.get_predictor(p)
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
            print ('predictor:', p)
            print ('alleles:',a)
            print ('length:',length)
            print ('cpus:', self.cpus)
            if 'iedb' in p:
                if iedb_checks(method) == False:
                    continue

            P.predictProteins(self.sequences, length=length, alleles=a, names=self.names,
                              path=savepath, overwrite=self.overwrite, verbose=self.verbose,
                              method=method, cpus=self.cpus)
            #load results into predictor
            P.load(path=savepath)
            if P.data is None:
                print ('no results were found, did predictor run?')
                return
            preds.append(P)
            #print (preds)
            cutoff = self.cutoff
            cutoff_method = self.cutoff_method
            n = self.n
            print ('-----------------------------')

        self.preds = preds
        self.analysis()
        if self.plots == True:
            self.plot_results()
        return

    def analysis(self):
        """Do analysis of predictions"""

        preds = self.preds
        for P in self.preds:
            print (P)
            p = P.name
            cutoff = self.cutoff
            cutoff_method = self.cutoff_method
            n = self.n
            b = P.getBinders(cutoff=cutoff, cutoff_method=cutoff_method)
            print ('%s/%s binders' %(len(b), len(P.data)))
            if len(b) == 0:
                print ('no binders found, check your cutoff value')
                return
            pb = P.promiscuousBinders(binders=b, n=int(n), cutoff=cutoff, cutoff_method=cutoff_method)
            print ('found %s promiscuous binders at cutoff=%s, n=%s' %(len(pb),cutoff,n))
            pb.to_csv(os.path.join(self.path,'prom_binders_%s_%s.csv' %(p,n)))
            if self.verbose == True and len(pb)>0:
                print ('top promiscuous binders:')
                print (pb[:10])
            if self.genome_analysis == True:
                cl = analysis.find_clusters(pb, genome=self.sequences)
                cl.to_csv(os.path.join(self.path,'clusters_%s.csv' %p))
            print ('-----------------------------')
        return

    def plot_results(self):
        """Plot results of predictions"""

        preds = self.preds
        prots = self.names
        import pylab as plt
        height = 2*len(preds)
        path = os.path.join(self.path, 'plots')
        if not os.path.exists(path):
            os.mkdir(path)
        #print (self.preds)
        for prot in prots:
            #ax = plotting.plot_tracks(preds,name=prot,n=1,cutoff=self.cutoff,
            #                              figsize=(14,height),legend=True)
            #plotting.mpl_plot_regions(coords, ax, color='gray')
            ax = plotting.plot_bars(preds[0], prot, cutoff=20, chunks=2)
            plt.tight_layout()
            plt.savefig(os.path.join(path,prot), dpi=150)
            plt.close()
        print ('saved plots')
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

def iedb_checks(method):
    if check_iedbmhc1_path() == False:
        return False
    if check_iedb_method(method) == False:
        return False
    return True

def check_mhc1_length(l):
    if l<9 or l>13:
        print ('use MHCI n-mer lengths from 9-13')
        return False

def check_iedbmhc1_path():
    if not os.path.exists(base.iedbmhc1path):
        print ('IEDB MHC-I tools not found, check path')
        return False

def check_iedbmhc2_path():
    if not os.path.exists(base.iedbmhc2path):
        print ('IEDB MHC-II tools not found, check path')
        return False

def check_iedb_method(method):
    P1 = base.IEDBMHCIPredictor()
    P2 = base.IEDBMHCIIPredictor()
    m = P1.methods+P2.methods
    if method not in m:
        print ('%s not found. following available:' %method)
        print (m)
        return False
    return True

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

def test_run():
    """Test run for a HIV sample file"""

    path = os.path.dirname(os.path.abspath(__file__))
    options = config.baseoptions
    options['base']['sequence_file'] = os.path.join(path, 'testing','HIV-1.fa')
    options['base']['mhc2_alleles'] = 'human_common_mhc2'
    options['base']['path'] = 'hiv1_test'
    options['base']['verbose'] = True
    options = config.check_options(options)
    W = WorkFlow(options)
    st = W.setup()
    W.run()

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
                        default=False, help="Do test predictions")
    parser.add_option("-a", "--analysis", dest="analysis",
                        help="Analysis path", metavar="FILE")
    parser.add_option("-s", "--server", dest="server",  action="store_true",
                        default=False, help="Run web app")
    parser.add_option("-x", "--port", dest="port", default=8888,
                        help="Port for web app, default 8888")
    parser.add_option("-v", "--version", dest="version", action="store_true",
                        help="Get version")

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
    elif opts.server == True:
        #from epitopepredict.server import webapp
        #webapp.run(port=5000, debug=True)
        import epitopepredict.tornado_serve
        epitopepredict.tornado_serve.main(opts.port)
    elif opts.test == True:
        test_run()
        print ('these test predictions can be viewed in the web app')
    elif opts.version == True:
        from . import __version__
        print ('epitopepredict version %s' %__version__)
    else:
        print_help()

if __name__ == '__main__':
    main()
