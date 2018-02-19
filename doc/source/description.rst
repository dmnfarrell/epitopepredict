Introduction
============

epitopepredict provides a standardized programmatic interface for executing several MHC binding prediction methods.
The results from each method can then be processed and visualized in a consistent manner. The Tepitope module
implements the TEPITOPEPan method and requires no external program to run. netMHCIIpan must be downloaded
separately from the website and installed on your system. The process is quite simple. The same applies for
the IEDB tools. Both of these tools are free for academic use.

Command Line Interface
===========================

Installing the package provides the command epitopepredict in your path. This is a command line interface
to the library without the need for any Python coding. It provides pre-defined functionality with settings
specified in a text configuration file. Using this you can make MHC predictions with your chosen alleles and
predictors. If you are using the IEDB prediction tools they should be installed locally and you can specify
the path in the [iedbtools] section. Otherwise ignore those settings. Note that if settings are left out
generally defaults will be used so you can have a minimal file as in the examples.

Links
=====

http://openresearchsoftware.metajnl.com/articles/10.5334/jors.94/

http://dmnfarrell.github.io/epitopepredict

Installation
============

**Dependencies**

* numpy
* pandas
* matplotlib
* numexpr

**Optional dependencies**

* statsmodels
* seaborn (requires scipy)

This software has been tested on Linux.

Install on pip using:

`pip install epitopepredict`

Dependencies pandas and biopython should be installed automatically.

You can install as a snap which will usually be the latest version. Snaps are supported on many linux distros. Install snapd with your package manager and then:

`sudo snap install epitopepredict`