from setuptools import setup
import sys,os

setup(
    name = 'epitopepredict',
    version = '0.2.0',
    description = 'Python package for MHC binding prediction',
    url='https://github.com/dmnfarrell/epitopepredict',
    license='GPL v3',
    author = 'Damien Farrell',
    author_email = 'farrell.damien[at]gmail.com',
    packages = ['epitopepredict'],
    package_data={'epitopepredict': ['mhcdata/*.csv','tepitope/*.txt',
				     'tepitope/pssm/*']},
    install_requires=['numpy>=1.5',
                      'pandas>=0.17',
                      'biopython>=1.5',
                      'future'],
    entry_points = {
        'console_scripts': [
            'epitopepredict=epitopepredict.app:main']
            },
    classifiers = ['Operating System :: OS Independent',
            'Programming Language :: Python :: 2.7',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'License :: OSI Approved :: Apache Software License',
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research'],
    keywords = ['epitope','mhc','immunology','biology'],
)
