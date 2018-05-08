#!/usr/bin/env python

"""
    Command line script for neo epitope prediction
    Created March 2018
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, subprocess
import shutil
import pickle
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

        if check_imports() == False:
            return
        check_ensembl()
        pd.set_option('display.width', 120)
        base.iedbmhc1path = self.iedbmhc1_path
        base.iedbmhc2path = self.iedbmhc2_path

        self.vcf_files = self.vcf_files.split(',')
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

        if type(self.cutoffs) is int or type(self.cutoffs) is float:
            self.cutoffs = [self.cutoffs]
        else:
            self.cutoffs = [float(i) for i in self.cutoffs.split(',')]

        if os.path.exists(self.names):
            self.names = read_names(self.names)
        elif self.names == '':
            self.names=None
        else:
            self.names = self.names.split(',')
        if not os.path.exists(self.path) and self.path != '':
            os.mkdir(self.path)
        return True

    def get_file_labels(self, files):
        l=OrderedDict()
        for f in files:
            if not os.path.exists(f):
                print ('no such file %s' %f)
                continue
            label = os.path.splitext(os.path.basename(f))[0]
            l[label] = {'filename':f}
        return l

    def run(self):
        """Run workflow for multiple samples and prediction methods."""

        print ('running neoepitope predictions')
        path = self.path
        overwrite = self.overwrite
        files = self.vcf_files
        preds = self.predictors
        labels = self.get_file_labels(files)
        cutoffs = self.cutoffs
        if len(cutoffs) < len(preds) :
            cutoffs = [cutoffs[0] for p in preds]

        for f in labels:
            print (f)
            infile = labels[f]['filename']
            #file to save variants to, if present we can skip
            eff_csv = os.path.join(path, 'variant_effects_%s.csv' %f)
            eff_obj = os.path.join(path, 'variant_effects_%s.pickle' %f)
            if not os.path.exists(eff_obj) or overwrite == True:
                #get variant effects for each file and then iterate over predictors
                variants = load_variants(vcf_file=infile)
                labels[f]['variants'] = len(variants)
                print ('getting variant effects')
                effects = get_variant_effects(variants, self.verbose)
                #serialize variant effects
                effects_to_pickle(effects, eff_obj)
            else:
                #else reload from last run
                effects = pickle.load(open(eff_obj,'r'))
            #save as table
            eff_data = effects_to_dataframe(effects)
            eff_data['sample'] = f
            eff_data.to_csv(eff_csv)

            i=0
            for predictor in self.predictors:
                outfile = os.path.join(path, 'results_%s_%s.csv' %(f,predictor))
                if os.path.exists(outfile) and overwrite == False:
                    continue
                if predictor in base.mhc1_predictors:
                    alleles = self.mhc1_alleles
                    length = self.mhc1_length
                else:
                    alleles = self.mhc2_alleles
                    length = self.mhc2_length
                seqs = get_mutant_sequences(effects=effects, length=length, verbose=self.verbose)

                res = predict_variants(seqs, alleles=alleles, length=length,
                                 predictor=predictor, path=self.path, verbose=self.verbose, cpus=1)
                res['label'] = f
                res.to_csv(outfile, index=False)

                #gets promiscuous binders based on the cutoff
                P = base.get_predictor(predictor)
                P.data = res
                pb = P.promiscuous_binders(n=1, keep_columns=True, cutoff=cutoffs[i])
                #pb['label'] = f
                print (pb[:20])
                pb.to_csv(os.path.join(path, 'binders_%s_%s.csv' %(f,predictor)), index=False)
                i+=1
                #peps = self_similarity(res, proteome="human_proteome")

        #combine results for multiple samples
        pd.DataFrame(labels).T.to_csv(os.path.join(path, 'sample_labels.csv'))
        print ('finished, results saved to %s' %path)
        return

    def combine_results(self, labels):
        """Put peptides from multiple files in one table"""

        res=[]
        for i,r in labels:
            df=pd.read_csv('results_%s_tepitope.csv' %r.filename)
            res.append(df)
        res=pd.concat(res)
        pd.pivot_table(res, index=['peptide'], columns=['label'], values='score')
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
    """Get priority variant effects from a set of variants loaded with
    varcode. Omits silent and noncoding effects.
    Returns:
        varcode variant effect object
    """

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
    #if type(eff) is varcode.effects.effect_classes.FrameShift:
    #    print (eff, eff.shifted_sequence)
    if mut is None:
        return
    vloc = eff.aa_mutation_start_offset
    if vloc is None or len(v.coding_genes) == 0:
        return
    if verbose==True:
        print (eff.transcript_id, eff.short_description, eff.variant.ref)
        print (v.contig)
        #print eff.to_dict()['variant'].to_dict()
    return eff

def get_variant_effects(variants, verbose=False):
    """Get priority effects from a list of variants.
    Returns:
        list of varcode variant effect objects"""

    effects=[]
    for v in variants:
        eff = get_variant_effect(v, verbose=verbose)
        effects.append(eff)
    print ('%s effects' %len(effects))
    return effects

def effects_to_pickle(effects, filename):
    """serialize variant effects collections"""
    pickle.dump(effects, open(filename, "wb"), protocol=2)
    return

def effects_to_dataframe(effects):
    x=[]
    for eff in effects:
        if eff is None:
            continue
        d=OrderedDict()
        d['gene_name'] = eff.gene_name
        d['transcript_id'] = eff.transcript_id
        wt = eff.original_protein_sequence
        mut = eff.mutant_protein_sequence
        vloc = eff.aa_mutation_start_offset
        d['aa_change'] = eff.short_description
        d['mutation'] = eff.variant.short_description
        d['variant_class'] = get_variant_class(eff.variant)
        #d['wt_sequence'] = wt
        #d['mut_sequence'] = mut
        x.append(d)
    df = pd.DataFrame(x)
    df['chr'] = df.apply(lambda x: x.mutation.split(' ')[0],1)
    return df

def peptides_from_effect(eff, length=11, peptides=True):
    """Get mutated peptides from an effect object.
    Returns:
        dataframe with peptides and variant info
    """

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

def get_mutant_sequences(variants=None, effects=None, reference=None, peptides=True,
                         drop_duplicates=True, length=11, verbose=False):
    """
    Get mutant proteins or peptide fragments from vcf or maf file.
    Args:
        peptides:
    Returns:
        pandas dataframe with mutated peptide sequence and source information
    """

    res = []
    if variants is not None:
        effects = get_variant_effects(variants, verbose)
    if effects is None:
        print ('no variant information')
        return

    for eff in effects:
        peps = peptides_from_effect(eff, length=length, peptides=peptides)
        res.append(peps)
    res = pd.concat(res).reset_index(drop=True)
    if drop_duplicates == True:
        res = res.drop_duplicates('peptide')
    print ('%s sequences/peptides from %s effects' %(len(res),len(effects)))
    return res

def load_variants(vcf_file=None, maf_file=None, max_variants=None):
    import varcode
    if vcf_file is not None:
        variants = varcode.load_vcf(vcf_file, allow_extended_nucleotides=True, max_variants=max_variants)
        f=vcf_file
    elif maf_file is not None:
        variants = varcode.load_maf(maf_file)
        f=maf_file
    print ('%s variants read from %s' %(len(variants),f))
    return variants

def predict_variants(data, length=9,  predictor='tepitope', alleles=[],
                     path='', verbose=False, cpus=1, predictor_kwargs={}):
    """
    Predict binding scores for mutated peptides from supplied variants.
    Args:
        data: pandas dataframe with peptide sequences derived from get_mutant_sequences
    """

    if not os.path.exists(path):
        os.mkdir(path)
    df=data
    P = base.get_predictor(predictor)
    print ('predicting mhc binding for %s peptides with %s' %(len(df), P.name))
    peps = P.predict_peptides(df.peptide, alleles=alleles, cpus=cpus, **predictor_kwargs)
    peps = peps.drop(['pos','name'],1)

    res = df.merge(peps, on='peptide')
    #print (len(b), len(res))
    #x = self_similarity(pred, proteome="human_proteome")
    return res

def show_predictors():
    for p in base.predictors:
        print(p)

def check_imports():
    try:
        import varcode
    except:
        print ('varcode required. please run pip install varcode')
        return False
    return True

def fetch_ensembl_release(path=None, release='75'):
    """get pyensembl genome files"""

    from pyensembl import Genome,EnsemblRelease
    if path is not None:
        os.environ['PYENSEMBL_CACHE_DIR'] = path
    #this call should download the files
    genome = EnsemblRelease(release, species='human')
    print ('pyensembl genome files cached in %s' %genome.cache_directory_path)
    genome.download()
    genome.index()

def check_ensembl():
    """Check pyensembl ref genome cached. Needed for running in snap"""

    #check if running inside a snap package so we can download
    #the genome files for pyensembl
    if os.environ.has_key('SNAP_USER_COMMON'):
        print ('running inside snap')
        spath = os.environ['SNAP_USER_COMMON']
        print ('checking for ref human genome')
        fetch_ensembl_release()
    return

def print_help():
    print ("""use -h to get options""")

def plot_variant_summary(data):
    from bokeh.plotting import figure
    from bokeh.charts import Donut
    d = Donut(df, label=['abbr', 'medal'], values='medal_count',
          text_font_size='8pt', hover_text='medal_count')
    return d

def test_run():
    """Test run for sample vcf file"""

    print ('neoepitope workflow test')
    path = os.path.dirname(os.path.abspath(__file__))
    options = config.baseoptions
    options['base']['predictors'] = 'tepitope,mhcflurry'
    options['base']['mhc1_alleles'] = 'HLA-A*01:01,HLA-A*02:01,HLA-A*03:01'
    options['base']['mhc2_alleles'] = 'HLA-DRB1*01:01,HLA-DRB1*04:01,HLA-DRB1*08:01,HLA-DRB1*09:01,HLA-DRB1*11:01'
    options['base']['path'] = 'neo_test'
    options['base']['mhc2_length'] = 11
    #options['base']['verbose'] = True
    options['neopredict']['vcf_files'] = os.path.join(path, 'testing','input.vcf')
    options = config.check_options(options)
    #print (options)
    W = NeoEpitopeWorkFlow(options)
    st = W.setup()
    W.run()
