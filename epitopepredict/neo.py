#!/usr/bin/env python

"""
    Command line script for neo epitope prediction
    Created March 2018
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, subprocess
import time
import pickle
import numpy as np
import pandas as pd
pd.set_option('display.width', 150)
pd.set_option('max_colwidth', 80)
from collections import OrderedDict
from epitopepredict import base, config, analysis, sequtils, peptutils, tepitope

defaultpath = os.getcwd()
sim_matrix = tepitope.get_matrix('pmbec')
metrics = ['score', 'matched_score', 'binding_diff','perc_rank',
               'wt_similarity','self_similarity', 'virus_similarity',
               'anchor','hydro', 'net_charge']

class NeoEpitopeWorkFlow(object):
    """Class for implementing a neo epitope workflow."""

    def __init__(self, opts={}):
        for i in opts:
            self.__dict__[i] = opts[i]
        return

    def setup(self):
        """Setup main parameters"""

        if check_imports() == False:
            return
        #check_ensembl()
        pd.set_option('display.width', 120)
        base.iedbmhc1path = self.iedbmhc1_path
        base.iedbmhc2path = self.iedbmhc2_path

        self.vcf_files = self.vcf_files.split(',')
        f = self.vcf_files[0]
        fileext = os.path.splitext(f)[1]
        if fileext == '.txt' and os.path.exists(f):
            print ('found file list')
            self.vcf_files = read_names(f)

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

        self.names=None

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
        start = time.time()
        path = self.path
        overwrite = self.overwrite
        files = self.vcf_files
        preds = self.predictors
        labels = self.get_file_labels(files)
        cutoffs = self.cutoffs
        if len(cutoffs) < len(preds) :
            cutoffs = [cutoffs[0] for p in preds]

        for f in labels:
            print ('sample name: %s' %f)
            infile = labels[f]['filename']
            #file to save variants to, if present we can skip
            eff_csv = os.path.join(path, 'variant_effects_%s.csv' %f)
            eff_obj = os.path.join(path, 'variant_effects_%s.pickle' %f)
            if not os.path.exists(eff_obj) or overwrite == True:
                #get variant effects for each file and then iterate over predictors
                variants = load_variants(vcf_file=infile)
                labels[f]['variants'] = len(variants)
                print ('getting variant effects')
                effects = get_variants_effects(variants, self.verbose)
                #serialize variant effects
                effects_to_pickle(effects, eff_obj)
            else:
                #else reload from last run
                effects = pickle.load(open(eff_obj,'rb'))
            #save effects as table
            eff_data = effects_to_dataframe(effects)
            eff_data['sample'] = f
            eff_data.to_csv(eff_csv)

            #get mutated peptides
            seqs = get_mutant_sequences(effects=effects, length=self.mhc1_length, verbose=self.verbose)
            #get similarities
            df = get_closest_matches(seqs, self.verbose, cpus=self.cpus)

            i=0
            for predictor in self.predictors:
                outfile = os.path.join(path, 'results_%s_%s.csv' %(f,predictor))
                if os.path.exists(outfile) and overwrite == False:
                    continue
                if predictor in base.mhc1_predictors:
                    alleles = self.mhc1_alleles
                    #length = self.mhc1_length
                else:
                    alleles = self.mhc2_alleles
                    #length = self.mhc2_length

                res = predict_binding(df, alleles=alleles, predictor=predictor,
                                       verbose=self.verbose, cpus=self.cpus)
                res['label'] = f
                res.to_csv(outfile, index=False)

                #gets promiscuous binders based on the cutoff
                #P = base.get_predictor(predictor)
                #P.data = res
                #pb = P.promiscuous_binders(n=1, keep_columns=True, cutoff=cutoffs[i])
                #pb['label'] = f
                #print (pb[:20])
                #pb.to_csv(os.path.join(path, 'binders_%s_%s.csv' %(f,predictor)), index=False)
                i+=1

            #combine results if multiple predictors?
            #combine_results()

        #combine results for multiple samples?

        #save sample labels
        pd.DataFrame(labels).T.to_csv(os.path.join(path, 'sample_labels.csv'))
        print ('finished, results saved to %s' %path)
        end = round(time.time()-start,1)
        print ('took %s seconds' %end)
        return

    def combine_samples(self, labels):
        """Put peptides from multiple files in one table"""

        res=[]
        for i,r in labels:
            df=pd.read_csv('results_%s_tepitope.csv' %r.filename)
            res.append(df)
        res=pd.concat(res)
        pd.pivot_table(res, index=['peptide'], columns=['label'], values='score')
        return

def pbmec_score(seq1, seq2):
    """Score with PBMEC matrix"""
    x=0
    try:
        for i in seq1:
            for j in seq2:
                x+=sim_matrix[i][j]
    except:
        return -1
    return x

def get_alleles(f):
    """Get input alleles"""

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

def variants_from_csv(csv_file, sample_id=None, reference=None):
    """Variants from csv file.
    
    Args:
        csv_file: csv file with following column names-
            chromosome, position, reference_allele, alt_allele, gene_name, transcript_id, sample_id
        sample_id: if provided, select variants only for this id
        reference: ref genome used for variant calling
    """

    from pyensembl import ensembl_grch38
    import varcode
    from varcode import Variant
    df = pd.read_csv(csv_file)
    variants=[]
    if sample_id != None and 'sample_id' in df.columns:
        df = df[df.sample_id==sample_id]
        df = df.drop_duplicates(['POS','REF','ALT'])
    for i,r in list(df.iterrows()):
        #print i
        v = Variant(contig=r.CHROM, start=r.POS, ref=r.REF, alt=r.ALT, ensembl=ensembl_grch38)
        variants.append(v)
    varcl = varcode.variant_collection.VariantCollection(variants)
    return varcl

def dataframe_to_vcf(df, outfile):
    """Write a dataframe of variants to a simple vcf file. Dataframe requires
    the following columns: #CHROM','POS','ID','REF','ALT'
    """

    f = open(outfile, 'w')
    f.write('##fileformat=VCFv4.0\n')
    f.write('##reference=GRCh38.fa\n')
    f.write('##source=csv\n')
    f.write('##FILTER=<ID=PASS,Description="Passed all filters">\n')
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    df = df.rename(columns={'CHROM':'#CHROM'})
    df['ID']='.'
    df['FILTER']='PASS'
    df['QUAL'] = 60
    df['INFO'] = '.'
    df['FORMAT'] = 'GT'
    df['sample'] = '0/1'
    cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','sample']
    df[cols].to_csv(f, sep='\t',index=False)#,align='left')
    return

def get_variant_class(effect):

    import varcode
    v = effect.variant
    if v.is_deletion:
        return 'deletion'
    elif v.is_insertion:
        return 'insertion'
    elif v.is_snv:
        return 'snv'
    elif v.is_indel:
        return 'indel'
    elif type(effect) is varcode.effects.effect_classes.FrameShift:
        return 'frameshift'

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
        d['variant_class'] = get_variant_class(eff)
        #d['wt_sequence'] = wt
        #d['mut_sequence'] = mut
        x.append(d)
    df = pd.DataFrame(x)
    df['chr'] = df.apply(lambda x: x.mutation.split(' ')[0],1)
    return df

'''def filter_variant_effects(effects, verbose=False):
    """Get fitered list of effects.
       Omits silent and noncoding effects.
       Returns:
        list of varcode variant effect objects
    """

    effs = effects.drop_silent_and_noncoding()
    filt = []
    for eff in effs:
        v = eff.variant
        mut = eff.mutant_protein_sequence
        #if type(eff) is varcode.effects.effect_classes.FrameShift:
        #    print eff, eff.shifted_sequence
        if mut is None:
            continue
        vloc = eff.aa_mutation_start_offset
        if vloc is None or len(v.coding_genes) == 0:
            continue
        if verbose==True:
            print (v.gene_names, type(eff))
            print (eff.transcript_id, eff.short_description, eff.variant.ref)
        #print v.contig
        #print mut
        #print eff.to_dict()['variant'].to_dict()
        filt.append(eff)
    return filt'''


def get_variants_effects(variants, verbose=False, gene_expression_dict=None):
    """Get all effects from a list of variants.
    Returns:
        list of varcode variant effect objects"""

    from varcode.effects import Substitution, Insertion, Deletion
    effects = variants.effects()
    effects = effects.drop_silent_and_noncoding()
    #print (len(effects))

    effects = effects.filter_by_effect_priority(Substitution)
    if gene_expression_dict is not None:
         effects = effects.filter_by_gene_expression(gene_expression_dict)

    print ('%s effects from %s variants' %(len(effects),len(variants)))
    return effects

def peptides_from_effect(eff, length=11, peptides=True, verbose=False):
    """Get mutated peptides from a single effect object.
    Returns:
        dataframe with peptides and variant info
    """

    import varcode
    pad = length-1
    if eff==None:
        return
    gene = eff.gene_name
    varclass = get_variant_class(eff)
    orig = eff.original_protein_sequence
    mut = eff.mutant_protein_sequence
    vloc = eff.aa_mutation_start_offset
    if vloc is None or mut is None:
        return
    st = vloc-pad; end = vloc+pad+1
    if st<0: st=0

    #if frameshift changed sequence may be long
    if type(eff) is varcode.effects.effect_classes.FrameShift:
        #mutpep = mut[vloc:]
        mutpep = eff.shifted_sequence
        #print (mutpep)
        wt = None
        #print (type(eff), len(orig), len(mut), vloc, len(mutpep), mutpep, mut[vloc:], wt)
        #print (mut)
    else:
        mutpep = mut[st:end]
        if varclass == 'snv':
            wt = orig[st:end]
        else:
            wt = None
    #if verbose == True:
    #    print (type(eff), len(orig), len(mut), vloc, st, end, len(mutpep))
    if len(mutpep)<length:
        #if verbose == True:
        #    print ('peptide length too small')
        return
    if peptides is True:
        df = peptutils.get_fragments(seq=mutpep, length=length)
        df['pos'] = pd.Series(range(st,end))
        df['prot_length_ratio'] = len(mut)/float(len(orig))
        if wt != None:
            wdf = peptutils.get_fragments(seq=wt, length=length)
            df['wt'] = wdf.peptide
        else:
            df['wt'] = None
        #print (df)
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
    df['variant_class'] = varclass
    return df

def get_mutant_sequences(variants=None, effects=None, reference=None, peptides=True,
                         drop_duplicates=True, length=11, verbose=False):
    """
    Get mutant proteins or peptide fragments from vcf or maf file.
    Args:
        variants: varcode variant collection
        effects: non-synonmymous effects, alternative to variants
        peptides: get peptide fragments around mutation
    Returns:
        pandas dataframe with mutated peptide sequence and source information
    """

    res = []
    if variants is not None:
        effects = get_variants_effects(variants, verbose)
    if effects is None:
        print ('no variant information')
        return

    for eff in effects:
        #print (eff)
        peps = peptides_from_effect(eff, length=length, peptides=peptides, verbose=verbose)
        if peps is None:
            continue
        res.append(peps)
    res = pd.concat(res).reset_index(drop=True)

    #remove rows where mut same as wt peptide
    res = res[res.peptide!=res.wt]
    if drop_duplicates == True:
        res = res.drop_duplicates('peptide')
    print ('%s sequences/peptides from %s effects' %(len(res),len(effects)))
    return res

def load_variants(vcf_file=None, maf_file=None, max_variants=None):
    """Load variants from vcf file"""

    import varcode
    if vcf_file is not None:
        variants = varcode.load_vcf(vcf_file, allow_extended_nucleotides=True, max_variants=max_variants)
        f=vcf_file
    elif maf_file is not None:
        variants = varcode.load_maf(maf_file)
        f=maf_file
    print ('%s variants read from %s' %(len(variants),f))
    return variants

def get_closest_matches(df, verbose=False, cpus=1):
    """
    Find peptide similarity metrics
    """

    if verbose == True:
        print ('finding matches to self proteome')
    #find matches to self proteome, adds penalty score column to df
    #we should only blast non-duplicates....
    df = self_matches(df, cpus=cpus)
    if verbose == True:
        print ('finding matches to viral proteomes')
    df = virus_matches(df, cpus=cpus)
    #get similarity scores for wt and closest match to proteome
    matrix = 'pmbec'
    df['wt_similarity'] = df.apply(lambda x: wt_similarity(x, matrix=matrix),1)
    df['self_similarity'] = df.apply(lambda x: self_similarity(x, matrix=matrix),1)
    df['virus_similarity'] = df.apply(lambda x: virus_similarity(x, matrix=matrix),1)
    #get closest peptide in another column, either wt or nearest self
    df['closest'] = df.apply(get_closest_match, 1)
    df['length'] = df.peptide.str.len()
    #exclude exact matches to self?
    if verbose==True:
        print ('%s peptides with exact matches to self' %len(df[df.self_mismatches==0]))
    return df

def predict_binding(df, predictor='netmhcpan', alleles=[],
                     verbose=False, cpus=1, cutoff=.95, cutoff_method='default'):
    """
    Predict binding scores for mutated and wt peptides (if present) from supplied variants.

    Args:
        df: pandas dataframe with peptide sequences, requires at least 2 columns
            'peptide' - the mutant peptide
            'wt' - a corresponding wild type peptide
        this data could be generated from get_mutant_sequences or from an external program
        predictor: mhc binding prediction method
        alleles: list of alleles
    Returns:
        dataframe with mutant and wt binding scores for all alleles
    """

    P = base.get_predictor(predictor, scoring='ligand')
    print (P)
    print ('predicting mhc binding for %s peptides with %s' %(len(df), P.name))

    peps = list(df.peptide)
    res = P.predict_peptides(peps, alleles=alleles, cpus=cpus,
                             cutoff=cutoff, cutoff_method=cutoff_method, drop_columns=True)

    if res is None:
        print ('no binding predictions!')
        return

    #predict closest matching peptide affinity
    if verbose == True:
        print ('predicting wt peptides')
    wtpeps = list(df.closest)
    #print wild type peptides
    b_wt = P.predict_peptides(wtpeps, alleles=alleles, cpus=cpus,
                               cutoff=cutoff, cutoff_method=cutoff_method, drop_columns=True)

    #combine mutant and matching binding predictions
    res = combine_wt_scores(res, b_wt, P.scorekey)
    res = res.drop(['pos','name'],1)

    #combine binding results with main dataframe
    res = df.merge(res, on='peptide')
    res['binding_diff'] = res[P.scorekey]/res.matched_score

    #anchor position mutated in any 9-mers
    res['anchor'] = res.apply(anchor_mutated, 1)
    #hydrophobicity and net charge
    res = analysis.peptide_properties(res, 'peptide')
    res['length'] = res.peptide.str.len()

    #merge promiscuity measure into results
    #if len(pb) > 0:
    #    res = res.merge(pb[['peptide','alleles']], on='peptide',how='left')
    #else:
    #    res['alleles'] = 0
    #rename some columns
    res = res.rename(columns={'rank':'binding_rank','alleles':'promiscuity'})
    res = res.sort_values('binding_rank', ascending=True)
    return res

def score_peptides(df, rf=None):
    """Score peptides with a classifier. Returns a prediction probability."""

    if rf is None:
        from sklearn.externals import joblib
        rf = joblib.load(os.path.join(base.datadir,'rf_model.joblib'))
    X = df[metrics]
    X = X.fillna(X.mean())
    X = X.replace(np.inf,.1)
    sc = rf.predict_proba(X)[:,1]
    return sc

def combine_wt_scores(x, y, key):
    """Combine mutant peptide and matching wt/self binding scores from a
    set of predictions. Assumes both dataframes were run with the same alleles.
    Args:
        x,y: pandas dataframes with matching prediction results
        key:
    """

    x = x.sort_values(['pos','allele']).reset_index(drop=True)
    y = y.sort_values(['pos','allele']).reset_index(drop=True)
    x['matched_score'] = y[key]
    return x

def make_blastdb(url, name=None, filename=None, overwrite=False):
    """Download protein sequences and a make blast db. Uses datacache module."""

    import datacache
    cachedir = datacache.get_data_dir()
    blastdb = os.path.join(cachedir, name)
    if os.path.exists(blastdb+'.phr') and overwrite==False:
        #print ('blast files found')
        return blastdb

    filename = datacache.fetch_file(url, filename=filename, decompress=True, subdir=None)
    #print filename
    cmd = 'makeblastdb -dbtype prot -in %s -out %s' %(filename,blastdb)
    #print cmd
    tmp=subprocess.check_output(cmd, shell=True)
    return blastdb

def make_human_blastdb():
    """Human proteome blastdb"""

    url = 'ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz'
    filename = 'Homo_sapiens.GRCh38.pep.all.fa.gz'
    blastdb = make_blastdb(url, name='GRCh38', filename=filename)
    return blastdb

def make_virus_blastdb():
    """Human virus blastdb"""

    url = 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query=taxonomy:%22Viruses%20[10239]%22%20\
    keyword:%22Reference%20proteome%20[KW-1185]%22%20host:%22Homo%20sapiens%20(Human)%20[9606]%22&fil=&force=no&preview=true&format=fasta'
    filename = 'uniprot_human_virus_proteome.fa.gz'
    blastdb = make_blastdb(url, name='human_virus', filename=filename)
    return blastdb

def self_matches(df, **kwargs):

    blastdb = make_human_blastdb()
    x = find_matches(df, blastdb, **kwargs)
    x = x.rename(columns={'sseq':'self_match','mismatch':'self_mismatches'})
    return x

def virus_matches(df, **kwargs):

    blastdb = make_virus_blastdb()
    x = find_matches(df, blastdb, **kwargs)
    if 'sseq' in x.columns:
        x = x.rename(columns={'sseq':'virus_match','mismatch':'virus_mismatches'})
    else:
        x['virus_match'] = None
    return x

def find_matches(df, blastdb, cpus=4, verbose=False):
    """
    Get similarity measures for peptides to a self proteome. Does a
    local blast to the proteome and finds most similar matches. These
    can then be scored.
    Args:
        df: dataframe of peptides
        blastdb: path to protein blastdb
    Returns:
        dataframe with extra columns: 'sseq','mismatch'
    """
    if verbose == True:
        print ('blasting %s peptides' %len(df))
    length = df.peptide.str.len().max()

    def check_mm(x):
        #corrected mismatches for shorter hits
        if x.length<length:
            return length-x.length+x.mismatch
        else:
            return x.mismatch

    bl = sequtils.blast_sequences(blastdb, df.peptide, evalue=200000, cpus=cpus,
                                  ungapped=True, gapopen=10, gapextend=2, qcov_hsp_perc=100,
                                  comp_based_stats=0)

    if len(bl) == 0:
        if verbose == True:
            print ('no hits found!')
        return df
    if verbose == True:
        print ('%s hits' %len(bl))
    cols = ['qseqid','sseq','mismatch']

    #ignore any hits with gaps
    bl = bl[(bl.gapopen==0)]# & (bl.length>=length)]

    #take longest hit with lowest e-value for each query
    bl = bl.sort_values(['qseqid','length','evalue'],ascending=(True,False,True))
    bl = bl.groupby(['qseqid'],as_index=False).first()

    #correct mismatches to account for shorter hits
    bl['mismatch'] = bl.apply(check_mm, 1)
    bl = bl[cols]
    #merge results
    x = df.merge(bl,left_on='peptide',right_on='qseqid', how='left')
    x = x.sort_values(by='mismatch',ascending=True)
    x = x.drop(['qseqid'],1)
    #x['exact_match'] = x.mismatch.clip(0,1).fillna(1)
    return x

def wt_similarity(x, matrix='blosum62'):

    x1 = x.peptide
    x2 = x.wt
    #print(x1,x2)
    matrix = tepitope.get_matrix(matrix)
    return tepitope.similarity_score(matrix,x1,x2)

def self_similarity(x, matrix='blosum62'):

    if x.self_match is None:
        return
    x1 = x.peptide
    x2 = x.self_match
    matrix = tepitope.get_matrix(matrix)
    return tepitope.similarity_score(matrix,x1,x2)

def virus_similarity(x, matrix='blosum62'):

    if x.virus_match is None:
        return
    x1 = x.peptide
    x2 = x.virus_match
    matrix = tepitope.get_matrix(matrix)
    return tepitope.similarity_score(matrix,x1,x2)

def get_closest_match(x):
    """Create columns with closest matching peptide.
    If no wt peptide use self match. vector method"""
    if x.wt is None:
        return x.self_match
    else:
        return x.wt

def anchor_mutated(x):
    return peptutils.compare_anchor_positions(x.wt, x.peptide)

def summary_plots(df):
    """summary plots for testing results"""

    f,axs=plt.subplots(2,2,figsize=(10,10))
    axs=axs.flat
    g = df.groupby(['name']).size().sort_values(ascending=False)[:20]
    g.plot(kind='barh',ax=axs[0],color='gray')
    axs[0].set_title('peptide counts')
    df.variant_class.value_counts().plot(kind='pie',autopct='%.1f',ax=axs[1])
    axs[1].set_title('variant classes')
    df.self_mismatches.value_counts().sort_index().plot(kind='bar',ax=axs[2])
    axs[2].set_title('mismatches to self')
    #df.wt_similarity.hist(ax=axs[3])
    #df.plot('wt_similarity','self_similarity',kind='scatter',ax=axs[3])
    df.plot('score','matched_score',kind='scatter',ax=axs[3])
    return

def show_predictors():
    for p in base.predictors:
        print(p)

def check_imports():
    try:
        import varcode
    except Exception as e:
        print (e)
        print ('varcode required. please run pip install varcode')
        return False
    return True

def fetch_ensembl_release(path=None, release='75'):
    """Get pyensembl genome files"""

    from pyensembl import Genome,EnsemblRelease
    #this call should download the files
    genome = EnsemblRelease(release, species='human')
    genome.download(overwrite=False)
    genome.index(overwrite=False)
    genome.cache_directory_path = path
    print ('pyensembl genome files cached in %s' %genome.cache_directory_path)
    #run_pyensembl_install()
    return

def check_ensembl(release='75'):
    """Check pyensembl ref genome cached. Needed for running in snap"""

    #check if running inside a snap package so we can download
    #the genome files for pyensembl
    cache_dir=None
    if base.check_snap() is True:
        #home = os.path.join('/home', os.environ['USER'])
        home = os.environ['SNAP_USER_COMMON']
        cache_dir = os.path.join(home, '.cache')
        os.environ['PYENSEMBL_CACHE_DIR'] = cache_dir
    print ('checking for ref human genome')
    fetch_ensembl_release(cache_dir, release)
    return

def run_vep(vcf_file, out_format='vcf', assembly='GRCh38', cpus=4, path=None):
    """Run ensembl VEP on a vcf file for use with pvacseq.
    see https://www.ensembl.org/info/docs/tools/vep/script/index.html
    """

    fname = os.path.splitext(vcf_file)[0]
    out = fname+'.vep.%s' %out_format
    if path == None:
        path = '/local/ensembl-vep/'
    path = os.path.join(path,'./vep')
    cmd = '{p} --input_file {i} --pick --force_overwrite \
    --assembly {a} --fork {c} \
    --symbol --terms SO --output_file {o} \
    --plugin Downstream --plugin Wildtype \
    --cache --offline'.format(o=out,i=vcf_file,a=assembly,c=cpus,p=path)
    if out_format == 'vcf':
        cmd += ' --format vcf --vcf'
    print (cmd)
    tmp = subprocess.check_output(cmd, shell=True)
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
    options['base']['predictors'] = 'netmhcpan,tepitope'
    options['base']['mhc1_alleles'] = 'HLA-A*02:01'
    options['base']['path'] = 'neo_test'
    options['base']['overwrite'] = True
    #options['base']['mhc2_length'] = 11
    #options['base']['verbose'] = True
    options['base']['cpus'] = 2
    options['neopredict']['vcf_files'] = os.path.join(path, 'testing','input.vcf')
    options['neopredict']['release'] = '75'
    options = config.check_options(options)
    #print (options)
    W = NeoEpitopeWorkFlow(options)
    check_ensembl(release='75')
    st = W.setup()
    #check_ensembl()
    W.run()

def varcode_test():
    path = os.path.dirname(os.path.abspath(__file__))
    infile = os.path.join(path, 'testing','input.vcf')
    variants = load_variants(vcf_file=infile)
    get_variants_effects(variants)
    return
