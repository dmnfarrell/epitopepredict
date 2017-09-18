#!/usr/bin/env python

"""
    MHC prediction unit tests
    Created September 2015
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os
import pandas as pd
import unittest
from . import base, analysis, sequtils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class PredictorTests(unittest.TestCase):
    """Basic tests for predictor"""

    def setUp(self):
        genbankfile = 'testing/zaire-ebolavirus.gb'
        self.df = sequtils.genbank_to_dataframe(genbankfile, cds=True)
        self.testdir = 'testing'
        if not os.path.exists(self.testdir):
            os.mkdir(self.testdir)
        return

    def test_tepitope(self):
        """Tepitope test"""

        df = self.df
        P = base.get_predictor('tepitope')
        alleles = ["HLA-DRB1*0101", "HLA-DRB1*0305"]
        print (P)
        P.predictProteins(df, length=11, alleles=alleles,
                          path=self.testdir)
        P.getBinders(data=P.data)
        return

    def test_netmhciipan(self):
        """netMHCIIpan test"""

        #requires netmHCIIpan is installed
        df = self.df
        P = base.get_predictor('netmhciipan')
        alleles = ["HLA-DRB1*0101"]
        names = ['ZEBOVgp1']
        print (P)
        P.predictProteins(df, length=11, alleles=alleles, names=names,
                          path=self.testdir)
        P.getBinders(data=P.data)
        return

    '''def testIEDB(self):
        """IEDB MHCI test"""

        df = self.df
        P = base.getPredictor('iedbmhc1')
        print (P)
        alleles = ["HLA-A*02:02", "HLA-A*11:01",
                   "HLA-B*15:17", "HLA-B*51:01",
                   "HLA-C*04:01", "HLA-E*01:03"]
        P.predictProteins(df, length=11, alleles=alleles,
                           path=self.testdir)
        return'''

    def test_bcell(self):
        """IEDB BCell test"""

        df = self.df
        names = ['VP24']
        P = base.get_predictor('iedbbcell')
        P.iedbmethod='Chou-Fasman'
        P.predictProteins(df, names=names, path=self.testdir)
        return

    def test_fasta(self):
        """Test fasta predictions"""

        fastafile = 'testing/zaire-ebolavirus.faa'
        df = sequtils.fasta_to_dataframe(fastafile)
        alleles = ["HLA-DRB1*0101"]
        P = base.get_predictor('tepitope')
        P.predictProteins(df, length=11, alleles=alleles, path=self.testdir)
        return

    def test_load(self):
        """Test re-loading predictions"""

        infile = os.path.join(self.testdir, 'ZEBOVgp1.csv')
        P = base.get_predictor('iedbmhc1')
        P.load(infile)
        return

    def test_save(self):
        """Test saving"""

        return

    def test_analysis(self):
        """Test analysis methods"""

        #get binders for an entire set of saved results
        #b = analysis.getAllBinders(path, method='tepitope', n=3)
        #find clusters of binders in these results
        #cl = analysis.findClusters(b, method, dist=9, minsize=3)
        return

    def test_features(self):
        """Test genbank feature handling"""

        fastafile = 'testing/zaire-ebolavirus.faa'
        df = sequtils.fasta_to_dataframe(fastafile)
        name = 'ZEBOVgp1'
        sequtils.dataframe_to_fasta(df)
        sequtils.check_tags(df)
        return

    def quit(self):
        self.app.quit()

def run():
    unittest.main()

if __name__ == '__main__':
    unittest.main()
