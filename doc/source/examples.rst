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
    #current options: 'tepitope', 'netmhciipan', 'iedbmhc1', 'iedbmhc2', 'mhcflurry', 'bcell'
    PT = base.getPredictor('tepitope')


get sequence data::

    #get data in genbank format into a dataframe
    df = sequtils.genbank2Dataframe(genbankfile, cds=True)
    #get sequences from fasta file
    df = sequtils.fasta2Dataframe(fastafile)


run predictions::

    #run a single sequence
    seq = ep.testsequence
    label = 'myprot' #optional label for your sequence
    PT.predict(sequence=seq, allele='HLA-DRB1*01:01', length=11, name=label)

    #run predictions for multiple proteins
    #for 2 alleles and save results to savepath
    alleles = ["HLA-DRB1*01:01", "HLA-DRB1*03:05"]
    PT.predictProteins(df,length=11,alleles=alleles,save=True,path=savepath)


load previous results into a predictor::

    PT.load(path=path) #where path stores csv files for multiple proteins
    PT.load(filename=file) # where file is a csv formatted file of prediction results (can be 1 or more proteins)



