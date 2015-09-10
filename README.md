
### Background

**mhcpredict** provides a standardized programmatic interface for executing several MHC binding prediction methods. The results from each method can then be processed and visualized in a consistent manner. The Tepitope module implements the TEPITOPEPan method and requires no external program to run. netMHCIIpan must be downloaded separately from the website and installed on your system. The process is quite simple. The same applies for the IEDB tools. Both of these tools are free for academic use.

####Supported methods:

* TEPITOPEPan 
* NetMHCIIpan
* IEDB MHCI tools
* IEDB BCell tools

### Dependencies

* pandas
* biopython

### Installation

Clone the git repository into your Python path. Not yet available on pypi.
To use netMHCIIpan you need in install and added the path of the executable to your PATH. If you get the following error: `bash: /local/bin/netMHCIIpan: /bin/tcsh: bad interpreter:`, it means you are missing tcsh and should install it with your package manager.

### Methodology

Predictors for each method inherit from the Predictor class and all implement a predict method for scoring a single sequence. This may wrap methods from other modules and/or call command line predictors. For example the TepitopePredictor uses the mhcpredict.tepitope module. This method should return a Pandas DataFrame. The predictProteins method is used for multiple proteins contained in a dataframe of sequences in a standard format. This is created from a genbank or fasta file (see examples below). For large numbers of sequences predictProteins should be called with save=True so that the results are saved as each protein is completed to avoid memory issues, since many alleles might be called for each protein. Results are saved with one file per protein in msgpack format.

The results are of the following form and are returned sorted by the score column.

```
     peptide       core  pos  score      name         allele  rank
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
```

where name is the protein identifier from the input file (a locus tag for example) and a score column which will differ between methods. MHC-II methods can be run for varying lengths, with the core usually being the highest scoring in that peptide/n-mer (but not always).

### Example usage

```
#import
from mhcpredict import base, sequtils, analysis
#get list of predictors
print base.predictors 

#get data in genbank format into a dataframe
df = sequtils.genbank2Dataframe(genbankfile, cds=True)
#get data in fasta format
df = sequtils.fasta2Dataframe(fastafile)

#create tepitope predictor
P = base.getPredictor('tepitope')

#run prediction for 2 alleles and save results to savepath
alleles = ["HLA-DRB1*0101", "HLA-DRB1*0305"]
P.predictProteins(df,length=11,alleles=alleles,save=True,path=savepath)
#read previous results
res = pd.read_msgpack(file)
#set this data for the predictor
#assumes the data is for the right predictor, need to add checks...
P.data = res

#get binders from existing results (res is a DataFrame)
P.getBinders(data=res)
#get binders for an entire set of related proteins, i.e. in a genome
b = analysis.getAllBinders(path, method='tepitope', n=3)
#find clusters of binders in these results
cl = analysis.findClusters(b, method, dist=9, minsize=3)
```                        
