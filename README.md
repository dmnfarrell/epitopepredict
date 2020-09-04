<img align="right" src=https://raw.githubusercontent.com/dmnfarrell/epitopepredict/master/img/logo.png width=150px>

[![PyPI version shields.io](https://img.shields.io/pypi/v/epitopepredict.svg)](https://pypi.python.org/pypi/epitopepredict/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/epitopepredict/badge/?version=latest)](https://epitopepredict.readthedocs.io/en/latest/?badge=latest)

### Background

**epitopepredict** provides a standardized programmatic interface and command line tool for executing multiple epitope prediction methods. Currently this largely consists of interfaces to several MHC binding prediction, the results of which can then be processed and visualized in a consistent manner. The Tepitope module implements the TEPITOPEPan method is provided as a 'built in' method. The IEDB tools and netMHCIIpan and mhcflurry are also supported. All of these tools are free for academic use. This software runs on most linux systems. Users are recommended to use the snap package for convenience. This software is under active development particularly with a view to improve the command line and web tools.

Documentation is at http://epitopepredict.readthedocs.io

### Installation

current release:

`pip install epitopepredict`

from github repository:

`pip install -e git+https://github.com/dmnfarrell/epitopepredict.git#egg=epitopepredict`

or for snap package:

`snap install epitopepredict`
