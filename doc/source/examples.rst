Code Examples
=============

This page is for those using the Python API. For those wanting to use the command line application see the
Command line interface page. General usage of this package is to provide quick access to binding prediction methods
and perform analysis on the results. There are multiple potential applications.

Basics
------

imports::

    import epitopepredict as ep
    from epitopepredict import base, sequtils, analysis, plotting

create a Predictor object::

    #get list of predictors
    print base.predictors
    ['tepitope', 'netmhciipan', 'iedbmhc1', 'iedbmhc2', 'mhcflurry', 'mhcnuggets', 'iedbbcell']
    p = base.get_predictor('tepitope')

get sequence data::

    #get data in genbank format into a dataframe
    df = sequtils.genbank2Dataframe(genbankfile, cds=True)
    #get sequences from fasta file
    df = sequtils.fasta2Dataframe(fastafile)

run predictions for a protein sequence::

    seq = ep.testsequence
    label = 'myprot' #optional label for your sequence
    p = base.get_predictor('tepitope')
    p.predict(sequence=seq, allele='HLA-DRB1*01:01', length=11, name=label)

run predictions for multiple proteins::

    #run for 2 alleles and save results to savepath
    alleles = ["HLA-DRB1*01:01", "HLA-DRB1*03:05"]
    p = base.get_predictor('tepitope')
    p.predict_proteins(df, length=11, alleles=alleles, save=True, path=savepath)

run predictions for a list of peptides::

    from epitopepredict import peptutils
    seqs = peptutils.create_random_sequences(5000)
    p = ep.get_predictor('tepitope')
    x = p.predict_peptides(seqs, alleles=alleles)

run with multiple cpus::

    x = p.predict_peptides(seqs, alleles=alleles, cpus=4)

load previous results into a predictor::

    p.load(path=path) #where path stores csv files for multiple proteins
    p.load(filename=file) # where file is a csv formatted file of prediction results (can be 1 or more proteins)

Analysis
--------

get all the binders using the current data loaded into the predictor::

    #default is to use percentile cutoff per allele, returns a dataframe
    p.get_binders(cutoff=5)

get binders for only one protein by top median rank::

    p.get_binders(name=name, cutoff=20, cutoff_method='rank')

get all promiscuous binders, returns a dataframe::

    pb = p.promiscuous_binders(n=2, cutoff=5)
    #same using score cutoff
    pb = p.promiscuous_binders(n=2, cutoff_method='score', cutoff=500)

find clusters of binders in these results::

    cl = analysis.find_clusters(b, method, dist=9, minsize=3)
