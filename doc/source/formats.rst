File Formats
============

The command line interface currently accepts protein sequences as fasta or genbank files. Users will probably be familiar with these formats anyway if they are using them. A few useful tips are provided here:

* try to use sensible names in fasta file headers as they are used as identifiers for the results. When you get fasta files from some sources they can have very long headers like this::

    >lcl|NC_001802.1_prot_NP_057849.4_1 [gene=gag-pol] [locus_tag=HIV1gp1] [db_xref=GeneID:155348] [protein=Gag-Pol] [exception=ribosomal slippage] [protein_id=NP_057849.4] [location=join(336..1637,1637..4642)] [gbkey=CDS]


By default the text before the first space is used as the identifier for each protein so that should be unique. In this case it will be `lcl|NC_001802.1_prot_NP_057849.4_1`. You can also include an option called `fasta_header_sep` in the configuration file that will split the fasta name with another symbol as well, in this way you can shorten the names further, but they should still be unique.

* only CDS (coding sequence) features are used from genbank files, as obviously these are the ones with sequences

* make sure the /translation qualifier is present in the features of the genbank file. Some files might not have it and therefore no sequence is present. A typical genbank feature looks like this::

     CDS             360172..360507
                     /locus_tag="lmo0332"
                     /experiment="EXISTENCE:[PMID:19448609]"
                     /note="lmo0332"
                     /codon_start=1
                     /transl_table=11
                     /product="hypothetical protein"
                     /protein_id="NP_463862.1"
                     /db_xref="GeneID:987567"
                     /translation="MIYYICALYTFISALVSFGFSLDALLKSRKVNGDALINAKYAVS
                     RSLSLLIVALGLFIFKSDAFLVALSLVMIGAQLFDGIIGIKISTFKTVGPLLTAVGNV
                     IMLILFLTI"
