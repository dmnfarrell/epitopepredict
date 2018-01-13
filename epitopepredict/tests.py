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

path = os.path.dirname(os.path.abspath(__file__))
testdir = os.path.join(path, 'testing')

class PredictorTests(unittest.TestCase):
    """Basic tests for predictor"""

    def setUp(self):
        self.genbankfile = os.path.join(testdir, 'zaire-ebolavirus.gb')
        self.fastafile = os.path.join(testdir, 'zaire-ebolavirus.faa')
        self.df = sequtils.genbank_to_dataframe(self.genbankfile, cds=True)
        self.testdir = testdir
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

    def test_iedbmhc1(self):
        """IEDB MHCI test"""

        df = self.df
        P = base.get_predictor('iedbmhc1')
        base.iedbmhc1path = '/local/iedbmhc1'
        print (P)
        if not os.path.exists(base.iedbmhc1path):
            print ('IEDB MHC-I not found')
            return
        alleles = ["HLA-A*02:02", "HLA-A*01:01"]
        for m in P.methods:
            if m == 'comblib_sidney2008': continue
            print (m)
            P.predictProteins(df, length=8, alleles=alleles,
                                    method=m)
            b = P.getBinders(data=P.data, cutoff=5, cutoff_method='rank')
            print ('%s binders' %len(b))
        return

    def test_iedbmhc2(self):
        """IEDB MHCII test"""

        df = self.df
        P = base.get_predictor('iedbmhc2')
        base.iedbmhc2path = '/local/iedbmhc2'
        print ()
        print (P)
        if not os.path.exists(base.iedbmhc1path):
            print ('IEDB MHC-II not found')
            return
        alleles = ["HLA-DRB1*01:01"]
        for m in P.methods:
            print (m)
            P.predictProteins(df, length=8, alleles=alleles,
                                    method=m)
            b = P.getBinders(data=P.data, cutoff=5)
            #print (P.data.score)
            print ('%s binders' %len(b))
        return

    '''def test_bcell(self):
        """IEDB BCell test"""

        df = self.df
        names = ['VP24']
        P = base.get_predictor('iedbbcell')
        P.iedbmethod='Chou-Fasman'
        P.predictProteins(df, names=names, path=self.testdir)
        return'''

    def test_fasta(self):
        """Test fasta predictions"""

        df = sequtils.fasta_to_dataframe(self.fastafile)
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

        df = sequtils.fasta_to_dataframe(self.fastafile)
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
