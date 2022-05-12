Code Examples
=============

This page is for those using the Python API. For those wanting to use the command line application see the
Command line interface page. General usage of this package is to provide convenient access to binding prediction methods
and perform analysis on the results. There are multiple potential applications.

Methodology
-----------

MHC binding and other prediction methods are implemented by inheriting from a `Predictor` object. All such classes should at minimum override the `predict` method for scoring a single sequence. This may wrap methods from other python packages or call command line predictors. For example the `TepitopePredictor` uses the `epitopepredict.tepitope` module provided with this package.

The predict method should return a Pandas DataFrame. The `predict_sequences` method is used for multiple protein sequences contained in a dataframe of sequences in a standard format. This is created from a genbank or fasta file (see examples below). For large numbers of sequences `predict_sequences` should be called with save=True so that the results are saved as each protein is completed to avoid memory issues, since many alleles might be called for each protein. Results are saved with one file per protein/sequence in csv format.

The results are of the following form and are returned sorted by the score column::

        peptide       core      pos  score      name         allele  rank
   198  VIFRLMRTNFL  FRLMRTNFL  198    3.4  ZEBOVgp1  HLA-DRB1*0101     1
   199  IFRLMRTNFLI  FRLMRTNFL  199    3.4  ZEBOVgp1  HLA-DRB1*0101     1
   200  FRLMRTNFLIK  FRLMRTNFL  200    3.4  ZEBOVgp1  HLA-DRB1*0101     1
   709  NRFVTLDGQQF  FVTLDGQQF  709    2.5  ZEBOVgp1  HLA-DRB1*0101     4
   710  RFVTLDGQQFY  FVTLDGQQF  710    2.5  ZEBOVgp1  HLA-DRB1*0101     4
   711  FVTLDGQQFYW  FVTLDGQQF  711    2.5  ZEBOVgp1  HLA-DRB1*0101     4
   70   DSFLLMLCLHH  FLLMLCLHH   70    2.0  ZEBOVgp1  HLA-DRB1*0101     7
   71   SFLLMLCLHHA  FLLMLCLHH   71    2.0  ZEBOVgp1  HLA-DRB1*0101     7
   72   FLLMLCLHHAY  FLLMLCLHH   72    2.0  ZEBOVgp1  HLA-DRB1*0101     7
   32   QGIVRQRVIPV  IVRQRVIPV   32    1.7  ZEBOVgp1  HLA-DRB1*0101    10

where name is the protein identifier from the input file (a locus tag for example) and a score column which will differ between methods. MHC-II methods can be run for varying lengths, with the core usually being the highest scoring in that peptide/n-mer (but not always).

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
    df = sequtils.genbank_to_dataframe(genbankfile, cds=True)
    #get sequences from fasta file
    df = sequtils.fasta_to_dataframe(fastafile)

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

run with multiple threads::

    x = p.predict_peptides(seqs, alleles=alleles, threads=4)

load previous results into a predictor::

    p.load(path=path) #where path stores csv files for multiple proteins
    p.load(filename=file) # where file is a csv formatted file of prediction results (can be 1 or more proteins)

Analysis
--------

get all the binders using the current data loaded into the predictor::

    #default is to use percentile cutoff per allele, returns a dataframe
    p.get_binders(cutoff=.95)

get binders for only one protein by top median rank::

    p.get_binders(name=name, cutoff=10, cutoff_method='rank')

get all promiscuous binders, returns a dataframe::

    pb = p.promiscuous_binders(n=2, cutoff=.95)
    #same using score cutoff
    pb = p.promiscuous_binders(n=2, cutoff_method='score', cutoff=500)

find clusters of binders in these results::

    cl = analysis.find_clusters(b, method, dist=9, minsize=3)
