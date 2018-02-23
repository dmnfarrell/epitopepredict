#!/usr/bin/env python

"""
    Command line script for neo epitope prediction
    Created March 2018
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os
import shutil
import pandas as pd
from collections import OrderedDict
from epitopepredict import base, config, analysis, sequtils, peptutils

defaultpath = os.getcwd()

class NeoEpitopeWorkFlow(object):
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

        peps = predict_variants(vcf_file=self.vcf_file, alleles=self.mhc2_alleles[:4], length=11,
                                predictor=self.predictors[0],
                                path=self.path, verbose=self.verbose, cpus=1)

        #peps = self_similarity(peps, proteome="human_proteome")
        return

def get_variant_class(v):
    if v.is_deletion:
        return 'deletion'
    elif v.is_insertion:
        return 'insertion'
    elif v.is_snv:
        return 'snv'
    elif v.is_indel:
        return 'indel'

def get_variant_effect(variant, verbose=False):

    import varcode
    v = variant
    if verbose==True:
        print (v, v.gene_names)
    effs = v.effects()
    effs = effs.drop_silent_and_noncoding()

    if len(effs)>0:
        eff = effs.top_priority_effect()
        if verbose==True:
            print (v, v.gene_names)
            print (get_variant_class(variant))
            print (type(eff))
            print (effs)
            print()
    else:
        return
    mut = eff.mutant_protein_sequence
    if type(eff) is varcode.effects.effect_classes.FrameShift:
        print (eff, eff.shifted_sequence)
    if mut is None:
        return
    #print mut
    vloc = eff.aa_mutation_start_offset
    if vloc is None or len(v.coding_genes) == 0:
        return
    if verbose==True:
        print (eff.transcript_id, eff.short_description, eff.variant.ref)
        print (v.contig)
        #print eff.to_dict()['variant'].to_dict()
    return eff

def get_mutated_peptides(eff, length=11, peptides=True):

    import varcode
    pad = length-1
    if eff==None:
        return
    gene = eff.gene_name
    orig = eff.original_protein_sequence
    mut = eff.mutant_protein_sequence
    vloc = eff.aa_mutation_start_offset
    st = vloc-pad; end = vloc+pad+1
    #print type(eff)
    if type(eff) is varcode.effects.effect_classes.FrameShift:
        #mutpep = mut[vloc:]
        mutpep = eff.shifted_sequence
        wt = ''
    else:
        wt = orig[st:end]
        mutpep = mut[st:end]
    if len(mutpep)<length:
        return
    if peptides is True:
        df = peptutils.get_fragments(seq=mutpep, length=length)
        df['pos'] = pd.Series(range(st,end))
    else:
        #just return the mutated protein
        df = pd.DataFrame.from_dict([{'wt_sequence':orig,'mutant_sequence': mut}])
        df['pos'] = vloc
    #print gene,st,end,mutpep
    #print df
    df['name'] = gene
    df['transcript_id'] = eff.transcript_id
    #df['transcript_name'] = eff.transcript_name
    df['aa_change'] = eff.short_description
    df['mutation'] = eff.variant.short_description
    df['variant_class'] = get_variant_class(eff.variant)
    return df

def get_mutant_sequences(vcf_file=None, maf_file=None, reference=None, peptides=True,
                         length=11, verbose=False, max_variants=None):
    """
    Get mutant proteins or peptide fragments from vcf or maf file.
    Args:
        peptides:
    Returns:
        pandas dataframe with mutated peptide sequence and source information
    """

    res = []

    import varcode
    if vcf_file is not None:
        variants = varcode.load_vcf(vcf_file, allow_extended_nucleotides=True, max_variants=max_variants)
    elif maf_file is not None:
        variants = varcode.load_maf(maf_file)
    print ('%s variants read' %len(variants))
    ecount = 0
    for v in variants:
        eff = get_variant_effect(v, verbose=verbose)
        peps = get_mutated_peptides(eff, length=length, peptides=peptides)
        res.append(peps)

    #res = pd.DataFrame(peptides, columns=['protein','peptide'])
    res = pd.concat(res).reset_index(drop=True)
    print ('%s sequences/peptides from %s variants' %(len(res),len(variants)))
    print (res.variant_class.value_counts())
    return res

def predict_variants(vcf_file, length=9, alleles=[], path='', verbose=False,
                     predictor='smm', max_variants=None, cpus=1,
                     predictor_kwargs={}):
    """predict variants pipeline"""

    if not os.path.exists(path):
        os.mkdir(path)
    seqs = get_mutant_sequences(vcf_file, length=length, verbose=verbose, max_variants=max_variants)
    seqs.to_csv(os.path.join(path,'variant_sequences.csv'))
    P = base.get_predictor(predictor)
    print ('predicting mhc binding for %s peptides with %s' %(len(seqs), P.name))
    pred = P.predict_peptides(seqs.peptide, alleles=alleles, cpus=cpus, **predictor_kwargs)
    pred = pred.drop(['pos','name'],1)
    res = seqs.merge(pred,on='peptide')
    print (len(pred), len(res))
    #x = self_similarity(pred, proteome="human_proteome")
    res.to_csv(os.path.join(path, 'results_%s.csv' %predictor))
    return res

def show_predictors():
    for p in base.predictors:
        print(p)

    def analysis(self):
        """Do analysis"""

        return

def check_imports():
    try:
        import varcode
    except:
        print ('please run pip install varcode')

def print_help():
    print ("""use -h to get options""")

def test_run():
    """Test run for sample vcf file"""

    print ('running test')
    path = os.path.dirname(os.path.abspath(__file__))
    options = config.baseoptions
    options['base']['mhc2_alleles'] = 'human_common_mhc2'
    options['base']['path'] = 'neo_test'
    #options['base']['mhc2_length'] = 11
    #options['base']['verbose'] = True
    options['neopredict']['vcf_file'] = os.path.join(path, 'testing','input.vcf')
    options = config.check_options(options)
    #print (options)
    W = NeoEpitopeWorkFlow(options)
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
                        default=False, help="Run the pipeline")
    parser.add_option("-t", "--test", dest="test",  action="store_true",
                        default=False, help="Do test predictions")
    parser.add_option("-v", "--version", dest="version", action="store_true",
                        help="Get version")

    opts, remainder = parser.parse_args()

    if opts.config != None:
        cp = config.parse_config(opts.config)
        options = config.get_options(cp)
        options = config.check_options(options)
    else:
        print ('no config file')
    if opts.test == True:
        test_run()
    elif opts.run == True:
        W = WorkFlow(options)
        st = W.setup()
        if st == True:
            W.run()
    elif opts.version == True:
        from . import __version__
        print ('epitopepredict version %s' %__version__)
    else:
        print_help()

if __name__ == '__main__':
    main()
