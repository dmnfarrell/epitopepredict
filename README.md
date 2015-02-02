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

### Installation

Clone the git repository into your Python path or run `pip install mhcpredict`.

### Usage
```
#import
from mhcpredict import base,genome,analysis
#get data
df = Genome.genbank2Dataframe(genbankfile, cds=True)
#create class
P = Base.getPredictor('tepitope')
#run prediction for 2 alleles and save results to savepath
alleles = ["HLA-DRB1*0101", "HLA-DRB1*0305"]
P.predictProteins(df,length=11,alleles=alleles,save=True,path=savepath)
#use previous results
df = pd.read_msgpack(file)
#get binders
P.getBinders(data=df)
#get binders for an entire set of saved results
b = analysis.getAllBinders(path, method='tepitope', n=3)
#find clusters of binders in these results
cl = analysis.findClusters(b, method, dist=9, minsize=3)
```                        
