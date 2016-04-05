
### Background

**epitopepredict** provides a standardized programmatic interface for executing several epitope prediction methods. Currently this largely consists of interfaces to several MHC binding prediction, the results of which can then be processed and visualized in a consistent manner. Many MHC binding predictors have been developed but usually not in an open source manner. The Tepitope module implements the TEPITOPEPan method is provided as a 'built in' method. netMHCIIpan must be downloaded separately from the website and installed on your system. The process is quite simple. The same applies for the IEDB tools. Both of these tools are free for academic use. It is hoped that other epitope predictors will be integrated.

####Supported methods:

* TEPITOPEPan 
* NetMHCIIpan
* IEDB MHCI tools
* IEDB BCell tools

### Dependencies

* pandas
* biopython

### Installation

`pip install epitopepredict`

To use netMHCIIpan you need in install and added the path of the executable to your PATH. If you get the following error: `bash: /local/bin/netMHCIIpan: /bin/tcsh: bad interpreter:`, it means you are missing tcsh and should install it with your package manager.
