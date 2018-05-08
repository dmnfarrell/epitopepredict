Command Line Interface
======================

Installing the package provides the command `epitopepredict` in your path. This is a command line interface
to the library without the need for any Python coding. It provides pre-defined functionality with settings
specified in a text configuration file. Using this you can make MHC predictions with your chosen alleles and
predictors. If you are using the IEDB prediction tools they should be installed locally and you can specify
the path in the [iedbtools] section. Otherwise ignore those settings. Note that if settings are left out
generally defaults will be used so you can have a minimal file as in the examples.

Usage
-----

Usage largely involves setting up the config file and having your input files prepared.
Running the command `epitopepredict -c <yourfilename>.conf` will create a new config file for you to work from if it doesn't exist.
Just edit this with a text editor and then to execute::

    epitopepredict -c <yourfilename>.conf -r

You can also test the pipeline after installing by running::

    epitopepredict -t

This will generate predictions using a set of sample HIV-1 sequences and save the results to a folder called hiv1_test which you can open in the web app to view (see below). This should work 'out of the box' as it only uses the built in prediction algorithm, tepitope.

Configuration file settings
---------------------------

The advantage of configuration files is in avoiding long commands that have to be remembered or are prone to mistakes. Also the config files can be kept to recall what setting we used or to copy them for another set of files. The current options available in the file are shown below::

    [base]
    predictors = tepitope
    mhc2_alleles = HLA-DRB1*01:01,HLA-DRB1*04:01
    mhc1_alleles = HLA-A*01:01
    mhc1_length = 11
    mhc2_length = 15
    n = 2
    cutoff_method = default
    cutoffs = .95
    sequence_file =
    path = results
    overwrite = no
    verbose = no
    names =
    cpus = 1

    [iedbtools]
    iedbmhc1_path =
    iedbmhc2_path =
    iedb_mhc1_method = IEDB_recommended
    iedb_mhc2_method = IEDB_recommended

Settings explained
------------------

+------------------+-----------------------------+------------------------------------------------------------------------------+
| name             | example value               | meaning                                                                      |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| predictors       | tepitope                    | name of predictor: e.g. tepitope, iedbmhc1, netmhciipan, mhcflurry           |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| mhc1_alleles     | HLA-A*01:01,HLA-A*03:01     | list of MHC-I alleles or preset name                                         |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| mhc2_alleles     | HLA-DRB1*0101,HLA-DRB1*0103 | list of MHC-II alleles or preset name                                        |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| mhc1_length      | 11                          | length of n-mers for MHC-I prediction                                        |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| mhc2_length      | 11                          | length of n-mers for MHC-II prediction                                       |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| n                | 3                           | minimum number of alleles for counting as promiscuous binders                |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| cutoff_method    | score                       | cutoff method: default, score or rank used for getting binders (see below)   |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| cutoffs          | .95                         | percentile/score/rank cutoff for counting binders                            |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| sequence_file    | zaire-ebolavirus.gb         | set of protein sequences in genbank or fasta format                          |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| peptide_file     | peptides.txt                | set of peptides in a plain text file, one per row                            |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| path             | results                     | folder to save results to, can be empty for current folder                   |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| overwrite        | no                          | overwrite the previous results                                               |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| names            | Rv0011c,Rv0019c             | subset of protein/sequence names to predict from your input file, optional   |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| verbose          | no                          | displays more information while running                                      |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| cpus             | 1                           | number of processors to use, use 0 for all available                         |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| iedbmhc1_path    |                             | folder where the IEDB MHC-I tools are installed, not required unless used    |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| iedbmhc2_path    |                             | folder where the IEDB MHC-II tools are installed, not required unless used   |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| iedb_mhc1_method | IEDB_recommended            | predictor to use within the IEDB MHC-I tools (see below)                     |
+------------------+-----------------------------+------------------------------------------------------------------------------+
| iedb_mhc2_method | IEDB_recommended            | predictor to use within the IEDB MHC-II tools (see below)                    |
+------------------+-----------------------------+------------------------------------------------------------------------------+

Cutoff methods
--------------

Methods for achieving an appropriate cutoff for considering a peptide to be a binder are somewhat arbitrary. They vary with the application. There are three methods provided to select binders:

* **default** - allele specific global cutoffs, this uses a percentile cutoff to select peptides using pre-calculated quantile scores for each allele. This may avoid an issue where certain alleles will dominate if using a single score cutoff. Though there is limited evidence to suggest this is more appropriate. Typical value would be .95 i.e. top 95% in each allele.
* **rank** - Select top ranking peptides in each sequence above the cutoff. This would be useful for small numbers of sequence but for a lot of proteins might produce too many false positives.
* **score** - Use a single score cutoff for all peptides/alleles. This is probably the standard method. Typical binding predictors produce an affinity score and a cutoff of 500 is used. However this might also produce a lot of false positives.

Binding promiscuity
-------------------

Promiscuous binders are those above the cutoffs in more than n alleles. The rationale for this is that a peptide is more likely to be immunogenic in your target population if it is a binder in multiple alleles. This may not be the case in reality of course. By default the command line tool will calculate the promiscuous binders to give you a unique list of peptides and include the number of alleles in which it is a binder. The table is ranked by this value and the maximum score over the alleles tested.

Preset allele lists
-------------------

For convenience there are some lists of common alleles that you can use without having to type allele names into the config file. These have been taken from various sources and are only a rough guide. Use `epitopepredict -p` to see the available presets. The format of allele names is discussed on the MHC Allele Nomenclature page.

The current selection is:

+---------------------+--------------------------------------------------------+
| name                | description                                            |
+---------------------+--------------------------------------------------------+
| mhc1_supertypes     | 6 MHC-I supertypes                                     |
+---------------------+--------------------------------------------------------+
| mhc2_supertypes     | 7 MHC-II supertypes                                    |
+---------------------+--------------------------------------------------------+
| us_caucasion_mhc1   | 30 most common US caucasion MHC-I                      |
+---------------------+--------------------------------------------------------+
| us_african_mhc1     | 30 most common US african MHC-I                        |
+---------------------+--------------------------------------------------------+
| human_common_mhc2   | 11 most prevalent HLA-DR alleles worldwide             |
+---------------------+--------------------------------------------------------+
| broad_coverage_mhc1 | 26 alleles providing broad coverage                    |
+---------------------+--------------------------------------------------------+
| bovine_like_mhc2    | 8 HLA-DR alleles chosen to approximate bovine response |
+---------------------+--------------------------------------------------------+

IEDB tool methods
-----------------

The IEDB combines multiple prediction methods into its tools. Generally it's recommended to use their consensus methods but individual methods may be preferred. You can specify these using the iedb_mhc*_method options. Remember they do not all support all alleles. See Installing IEDB MHC tools.

MHC-I::

    ann
    comblib_sidney2008
    consensus
    IEDB_recommended
    netmhcpan
    smm
    smmpmbec

MHC-II::

    comblib
    consensus3
    IEDB_recommended
    NetMHCIIpan
    nn_align
    smm_align
    sturniolo

Examples
--------

**MHC-II binding predictions for preset alleles of proteins in a genbank file**

Using preset allele lists saves you the trouble of writing the alleles out. You can get the built-in presets by using -p at the command line. If you provide MHC-I alleles for a class II predictor like tepitope the program will give an error. More cpus means speed improvements::

    [base]
    predictors = tepitope
    mhc2_alleles = human_common_mhc2
    n = 2
    cutoffs = .95
    sequence_file = zaire-ebolavirus.gb
    path = results
    cpus = 2

**A small set of peptides**

Say we want to predict for small list of peptides with multiple prediction methods and select the top 10 ranking in at least 3 alleles. Here input.txt is just simple text file with all the individual peptides. They should be of an appropriate length::

    [base]
    predictors = tepitope,mhcflurry
    mhc1_alleles = human_common_mhc2
    mhc2_alleles = human_common_mhc2
    cutoff_method = rank
    cutoffs = 10
    n=3
    path = results
    peptide_file = input.txt

**Strict cutoffs**

For selection you can use very strict score cutoff level or high global percentile. In this example we use a score cutoff so must provide a cutoff value for each method::

    [base]
    predictors = tepitope,netmhciipan
    mhc1_alleles = human_common_mhc2
    cutoff_method = score
    cutoffs = 6,50
    n=3
    path = results
    peptide_file = input.txt

Outputs
-------

In each results folder you will find a sub-folder for each method. This has csv files with the predictions for each sequence, if using multiple protein sequences. This is the primary raw output. These folders can be re-used as input in the analysis section without re-running predictions and read by the web interface for presentation if needed. There are also files of the form final_method_n.csv which contain the promiscuous binders for each method.
