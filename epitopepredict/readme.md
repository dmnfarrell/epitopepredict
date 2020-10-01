<img src=https://raw.githubusercontent.com/dmnfarrell/epitopepredict/master/img/logo.png width=150px>

[![PyPI version shields.io](https://img.shields.io/pypi/v/epitopepredict.svg)](https://pypi.python.org/pypi/epitopepredict/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/epitopepredict/badge/?version=latest)](https://epitopepredict.readthedocs.io/en/latest/?badge=latest)

### Background

**epitopepredict** provides a standardized programmatic interface and command line tool for executing multiple epitope prediction methods. Currently this largely consists of interfaces to several MHC binding prediction, the results of which can then be processed and visualized in a consistent manner. There is a built-in method for MHC-class I prediction and the TEPITOPEPan method is provided as a 'built in' method for MHC-class II. The IEDB tools and netMHCpan, netMHCIIpan and MHCFlurry are also supported. Those tools are free for academic use but have to be installed separately. This software runs on most linux systems. 

Documentation is at http://epitopepredict.readthedocs.io

### Installation

`pip install epitopepredict`

or for latest version on github:

`pip install -e git+https://github.com/dmnfarrell/epitopepredict.git#egg=epitopepredict`
