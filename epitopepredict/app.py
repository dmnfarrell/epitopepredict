#!/usr/bin/env python

"""
    MHC prediction command line script
    Created March 2016
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os
import numpy as np
import shutil
import pandas as pd
from collections import OrderedDict
from epitopepredict import base, config, analysis, sequtils, plotting, neo

defaultpath = os.getcwd()

class WorkFlow(object):
    """Class for implementing a prediction workflow from a set of options"""
    def __init__(self, opts={}):
        for i in opts:
            self.__dict__[i] = opts[i]
        return

    def setup(self):
        """Setup main parameters"""

        if base.check_snap() == True:
            print ('running in a snap')
            cmd = base.set_netmhcpan_cmd()
            print ('netmhcpan cmd set to:', cmd)

        pd.set_option('display.width', 120)
        #override base.defaults entries if provided in conf
        set_defaults(self.__dict__)
        self.sequences = None
        self.peptides = None
        if self.sequence_file is not '':
            self.sequences = get_sequences(self.sequence_file, header_sep=self.fasta_header_sep)
            print ('input is %s protein sequences' %len(self.sequences))
        elif self.peptide_file is not '':
            self.peptides = read_names(self.peptide_file)
            print ('input is %s peptides' %len(self.peptides))
        if self.sequences is None and self.peptides is None:
            print ('no valid inputs found')
            return False
        self.mhc1_alleles = get_alleles(self.mhc1_alleles)
        self.mhc2_alleles = get_alleles(self.mhc2_alleles)
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
        if self.mhc2_alleles[0] in base.mhc2_presets:
            self.mhc2_alleles = base.get_preset_alleles(self.mhc2_alleles[0])

        if type(self.cutoffs) is int:
            self.cutoffs = [self.cutoffs]
        else:
            self.cutoffs = [float(i) for i in self.cutoffs.split(',')]

        if os.path.exists(self.names):
            self.names = read_names(self.names)
        elif self.names == '':
            self.names=None
        else:
            self.names = self.names.split(',')
        if self.names is not None:
            print ('selected sequences:', self.names)

        if not os.path.exists(self.path) and self.path != '':
            os.mkdir(self.path)

        #copy input seqs to path
        if self.sequences is not None:
            self.sequences.to_csv(os.path.join(self.path, 'input.csv'),index=False)
        return True

    def run(self):
        """Run prediction workflow"""

        preds = []
        for p in self.predictors:
            P = base.get_predictor(p)
            if P.check_install() == False:
                print ('%s not installed' %P.name)
                continue
            savepath = os.path.join(self.path, p)
            if self.overwrite == True and os.path.exists(savepath):
                shutil.rmtree(savepath)
            if p in base.mhc1_predictors:
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
            print ('alleles:',', '.join(a))
            print ('length:',length)
            print ('cpus:', self.cpus)
            if 'iedb' in p:
                if iedb_checks(method) == False:
                    continue

            if self.peptides is not None:
                P.predict_peptides(self.peptides, length=length, alleles=a,
                                path=self.path, overwrite=self.overwrite, verbose=self.verbose,
                                method=method, cpus=self.cpus, compression=self.compression)
                #load the results into the predictor
                P.load()
            else:
                #print (savepath)
                #data saved to disk to avoid large memory usage
                P.predict_proteins(self.sequences, length=length, alleles=a, names=self.names,
                                path=savepath, overwrite=self.overwrite, verbose=self.verbose,
                                method=method, cpus=self.cpus, compression=self.compression)
                #keep reference to path where results saved
                P.path = savepath
                #clear data as we will reload during analysis from disk
                P.data = None
            #if P.data is None:
            #    print ('no results were found, did predictor run?')
            #    return
            preds.append(P)
            n = self.n
            print ('-----------------------------')

        self.preds = preds
        self.analysis()
        #if self.plots == True:
        #    self.plot_results()
        print ('results saved in the folder %s' %self.path)
        return

    def analysis(self):
        """Do analysis of prediction results."""

        preds = self.preds
        cutoffs = self.cutoffs
        if len(cutoffs) < len(preds) :
            cutoffs = [cutoffs[0] for p in preds]
        cutoff_method = self.cutoff_method
        i=0
        prom_binders = {}
        print ('analysing results..')
        for P in self.preds:
            p = P.name
            cutoff = cutoffs[i]
            n = self.n
            print (P.path)
            if P.data is not None:
                b = P.get_binders(cutoff=cutoff, cutoff_method=cutoff_method)
            elif P.path is not None:
                b = P.get_binders(path=P.path, cutoff=cutoff, cutoff_method=cutoff_method)
            else:
                print ('empty results?')
                continue
            if b is None:
                continue
            print ('%s: %s binders found' %(P, len(b)))
            if len(b) == 0:
                print ('no binders found, check your cutoff value')
                return

            pb = P.promiscuous_binders(binders=b, n=n, cutoff=cutoff, cutoff_method=cutoff_method)
            print ('found %s promiscuous binders at cutoff=%s, n=%s' %(len(pb),cutoff,n))
            pb.to_csv(os.path.join(self.path,'final_%s_%s.csv' %(p,n)), float_format='%g')
            prom_binders[p] = pb

            if len(pb)>0:
                print ('top promiscuous binders:')
                print (pb[:10])
            else:
                continue
            if self.sequences is not None:
                x = analysis.create_nmers(pb, self.sequences, how='split', length=20)
                x = analysis.peptide_properties(x, 'nmer')
                x.to_csv(os.path.join(self.path,'final_%s_%s.csv' %(p,n)), float_format='%g')
                #do further analysis if using protein sequences
                cl = analysis.find_clusters(pb)
                if len(cl) > 0:
                    #make peptide lists
                    cl = analysis.create_nmers(cl, self.sequences, how='split', length=20)
                    cl = analysis.peptide_properties(cl, 'nmer')
                    cl.to_csv(os.path.join(self.path,'clusters_%s.csv' %p))
                #make summary table
                summary = self.get_summary(P, pb, self.sequences, clusters=cl)
                summary.to_csv(os.path.join(self.path,'summary_%s.csv' %p))
            print ('-----------------------------')
            i+=1

        if len(prom_binders) > 1:
            print ('finding combined list of binders')
            comb = self.combine(prom_binders)
            comb.to_csv(os.path.join(self.path,'combined.csv'), float_format='%g')
        return

    def get_summary(self, P, pb, seqs, clusters=None):
        """Get summary table for sequence based predictions.
        """

        rmcols = ['type','product','pseudo']
        try:
            seqs = seqs.drop(rmcols, 1)
        except:
            pass
        x = pb.groupby('name').agg({'peptide':[base.first,np.size],
                                    P.scorekey:base.first}).reset_index()
        x.columns = [col[1]+'_'+col[0] if col[1]!='' else col[0] for col in x.columns.values]
        x = x.rename(columns={'size_peptide':'binders','first_%s' %P.scorekey:'max_score',
                              'first_peptide':'top_peptide'})

        if clusters is not None and len(clusters)>0:
            clusters = clusters.groupby('name').size().rename('clusters').to_frame().reset_index()
            x = x.merge(clusters,on='name',how='left')
        if seqs is not None:
            if not 'length' in seqs.columns:
                seqs['length'] = seqs.translation.str.len()
            x = x.merge(seqs,left_on='name',right_on='locus_tag')
            x['binder_density'] = (x.binders/x.length).round(3)

        x = x.sort_values(by='binder_density',ascending=False)
        #print (x)
        return x

    def combine(self, data):
        """Combine peptide binders present in all methods."""

        cols = ['peptide','alleles','mean']
        keys = list(data.keys())
        n1 = keys[0]
        df = data[n1][['name']+cols]
        for n in keys[1:]:
            if data[n] is None:
                continue
            df2 = data[n][cols]
            df = df.merge(df2, on=['peptide'], how='inner', suffixes=['_'+n1,'_'+n])
        return df

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

def get_alleles(f):
    """Get input alleles as text file or list"""

    fileext = os.path.splitext(f)[1]
    if fileext == '.txt' and os.path.exists(f):
        items = read_names(f)
    else:
        items = f.split(',')
    return items

def read_names(filename):
    """read plain text file of items"""

    with open(filename) as f:
        p = f.readlines()
    p = [x.strip() for x in p]
    p = list(filter(None, p))
    return p

def get_sequences(filename, header_sep=None):
    """Determine file type and get sequences"""

    ext = os.path.splitext(filename)[1]
    if ext in ['.fa','.faa','.fasta']:
        seqs = sequtils.fasta_to_dataframe(filename, header_sep=header_sep)
        #print ('found fasta file')
    elif ext in ['.gb','.gbk','.genbank','.gbff']:
        seqs = sequtils.genbank_to_dataframe(filename, cds=True)
    return seqs

def set_defaults(d):
    """Override default paths if provided in conf"""
    for key in ['iedbmhc1_path','iedbmhc2_path']:
        if key in d and d[key] != '':
            base.defaults[key] = d[key]

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
    P = base.get_predictor('iedbmhc1')
    if not os.path.exists(P.iedbmhc1_path):
        print ('IEDB MHC-I tools not found, check path')
        return False

def check_iedbmhc2_path():
    P = base.get_predictor('iedbmhc2')
    if not os.path.exists(P.iedbmhc2_path):
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
    n=7
    for p in base.predictors:
        print (p)
        print ('-----------------------------')
        P = base.get_predictor(p)
        x = P.get_alleles()
        if type(x) is list:
            for i in range(0, len(x), n):
                l=x[i:i+n]
                print(', '.join(l))

        print ()
    return

def add_path():
    """Add home dir to path for accessing tools from a snap"""

    #toolspath = os.environ['SNAP_USER_COMMON']
    toolspath = os.path.join('/home', os.environ['USER'])
    binpath = os.path.join(toolspath, 'bin')
    if not os.path.exists(binpath):
        os.mkdir(binpath)
    print ('you should install external tools in %s' %toolspath)
    os.environ["PATH"] += os.pathsep + binpath
    #print (os.environ["PATH"])
    return

def check_installed():
    """Check which predictors can be used"""

    cl=base.get_predictor_classes()
    found=[]
    for i in cl:
        P=base.get_predictor(i)
        st = P.check_install()
        if st is True:
            found.append(i)
    return found

def test_run():
    """Test run for a sample file."""

    #installed = ','.join(check_installed())
    installed = 'tepitope,mhcflurry'
    path = os.path.dirname(os.path.abspath(__file__))
    options = config.baseoptions
    b=options['base']
    b['sequence_file'] = os.path.join(path, 'testing','zaire-ebolavirus.faa')
    b['mhc1_alleles'] = 'HLA-A*02:01,HLA-A*03:01,HLA-A*02:06'
    b['mhc2_alleles'] = 'human_common_mhc2'
    b['mhc1_length'] = 9
    b['predictors'] = installed
    b['path'] = 'zaire_test'
    b['verbose'] = True
    b['cutoff_method'] = 'score'
    b['cutoffs'] = '5,500'
    b['overwrite'] = True
    options = config.check_options(options)
    W = WorkFlow(options)
    st = W.setup()
    W.run()

def test_binary():

    cmd = base.set_netmhcpan_cmd()
    cmd = cmd + ' -a HLA-A02:01 -l 9 P170919-test.fa'
    print (cmd)
    import subprocess
    tmp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    import io
    tmp = tmp.decode()
    for i in tmp.split('\n'):
        print (i)
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
                        default=False, help="Do test predictions")
    parser.add_option("-a", "--analysis", dest="analysis",
                        help="Analysis path", metavar="FILE")
    parser.add_option("-n", "--neoepitope", dest="neoepitope", action="store_true",
                        default=False, help="Neo-epitope pipeline")
    parser.add_option("-e", "--ensembl", dest="ensembl", action="store_true",
                        default=False, help="Get ensembl files for a release")
    parser.add_option("-s", "--server", dest="server", action="store_true",
                        default=False, help="Run web app")
    parser.add_option("-x", "--port", dest="port", default=8000,
                        help="Port for web app, default 8000")
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
            config.write_config(conffile, defaults=config.baseoptions)
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
    elif opts.neoepitope == True:
        if opts.test == True:
            neo.test_run()
            #neo.varcode_test()
        else:
            print (options)
            release = options['ensembl_release']
            neo.check_ensembl(release)
            W = neo.NeoEpitopeWorkFlow(options)
            st = W.setup()
            if st == True:
                W.run()
    elif opts.server == True:
        #from epitopepredict.server import webapp
        #webapp.run(port=5000, debug=True)
        import epitopepredict.tornado_serve
        epitopepredict.tornado_serve.main(opts.port)
    elif opts.test == True:
        #test_binary()
        test_run()
    elif opts.version == True:
        from . import __version__
        print ('epitopepredict version %s' %__version__)
    else:
        print_help()

if __name__ == '__main__':
    main()
