mhcpredict
==========

MHC binding prediction tools scripts

### Background

This library provides a standardized programmatic interface for executing several MHC binding prediction methods. The results from each method can then be processed and visualized in a consistent manner. The Tepitope module implements the TEPITOPEPan method and requires no external program to run. netMHCIIpan must be downloaded separately from the website and installed on your system. The process is quite simple. The same applies for the IEDB tools. Both of these tools are free for academic use.

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

#get binders
P.getBinders(data=pred)
#get binders for an entire set of saved results
b = analysis.getAllBinders(path, method='tepitope', n=3)
#find clusters of binders in these results
cl = analysis.findClusters(b, method, dist=9, minsize=3)
```                        
