#!/usr/bin/env python

"""
    MHC prediction unit tests
    Created September 2015
    Copyright (C) Damien Farrell
"""

import sys, os
import pandas as pd
import unittest
import base, analysis, sequtils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class PredictorTests(unittest.TestCase):
    """Basic tests for predictor"""

    def setUp(self):
        genbankfile = 'testing/zaire-ebolavirus.gb'
        self.df = sequtils.genbank2Dataframe(genbankfile, cds=True)
        self.testdir = 'testing'
        if not os.path.exists(self.testdir):
            os.mkdir(self.testdir)
        return

    def testTepitope(self):
        """Tepitope test"""

        df = self.df
        P = base.getPredictor('tepitope')
        alleles = ["HLA-DRB1*0101", "HLA-DRB1*0305"]
        print P
        P.predictProteins(df, length=11, alleles=alleles,
                          save=True, path=self.testdir)
        P.getBinders(data=P.data)
        return

    def testnetMHCIIpan(self):
        """netMHCIIpan test"""

        #requires netmHCIIpan is installed
        df = self.df
        P = base.getPredictor('netmhciipan')
        alleles = ["HLA-DRB1*0101"]
        names = ['ZEBOVgp1']
        print P
        P.predictProteins(df, length=11, alleles=alleles, names=names,
                          save=True, path=self.testdir)
        P.getBinders(data=P.data)
        return

    '''def testIEDB(self):
        """IEDB MHCI test"""

        df = self.df
        P = base.getPredictor('iedbmhc1')
        print P
        alleles = ["HLA-A*02:02", "HLA-A*11:01",
                   "HLA-B*15:17", "HLA-B*51:01",
                   "HLA-C*04:01", "HLA-E*01:03"]
        P.predictProteins(df, length=11, alleles=alleles,
                          save=True, path=self.testdir)
        return'''

    def testBcell(self):
        """IEDB BCell test"""

        df = self.df
        names = ['VP24']
        P = base.getPredictor('bcell')
        P.iedbmethod='Chou-Fasman'
        P.predictProteins(df, names=names, save=True, path=self.testdir)
        return

    def testFasta(self):
        """Test fasta predictions"""

        fastafile = 'testing/zaire-ebolavirus.faa'
        df = sequtils.fasta2Dataframe(fastafile)
        alleles = ["HLA-DRB1*0101"]
        P = base.getPredictor('tepitope')
        P.predictProteins(df, length=11, alleles=alleles,
                          save=True, path=self.testdir)
        return

    def testLoad(self):
        """Test re-loading predictions"""

        infile = os.path.join(self.testdir, 'ZEBOVgp1.mpk')
        pred = pd.read_msgpack(infile)
        P = base.getPredictor('iedbmhc1')
        P.data = pred
        return

    def testAnalysis(self):
        """Test analysis methods"""

        #get binders for an entire set of saved results
        #b = analysis.getAllBinders(path, method='tepitope', n=3)
        #find clusters of binders in these results
        #cl = analysis.findClusters(b, method, dist=9, minsize=3)
        return

    def testFeatures(self):
        """Test genbank feature handling"""

        df = self.df
        name = 'ZEBOVgp1'
        #row = df[df.locus_tag==name]
        sequtils.dataframe2Fasta(df)
        sequtils.checkTags(df)
        #sequtils.fastaFormatfromFeature(feat)
        return

    def quit(self):
        self.app.quit()


if __name__ == '__main__':
    unittest.main()
