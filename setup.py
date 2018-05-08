from setuptools import setup
import sys,os

with open('epitopepredict/description.txt') as f:
    long_description = f.read()

setup(
    name = 'epitopepredict',
    version = '0.3.0',
    description = 'Python package for epitope prediction',
    long_description = long_description,
    url='https://github.com/dmnfarrell/epitopepredict',
    license='GPL v3',
    author = 'Damien Farrell',
    author_email = 'farrell.damien@gmail.com',
    packages = ['epitopepredict'],
    package_data={'epitopepredict': ['mhcdata/*.csv', 'presets/*.csv',
                  'tepitope/*.txt', 'tepitope/pssm/*', 'testing/*',
                  'templates/*.html','static/*',
                  'description.txt']
                 },
    install_requires=['numpy>=1.10',
                      'pandas>=0.20',
                      'biopython>=1.5',
                      'bokeh==0.12.14',
                      'wtforms>=2.1',
                      'wtforms_tornado',
                      'future'],
    entry_points = {
        'console_scripts': [
            'epitopepredict=epitopepredict.app:main']
            },
    classifiers = ['Operating System :: OS Independent',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.6',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords = ['epitope','mhc','immunology','biology'],
)
