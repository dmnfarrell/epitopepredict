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
from . import base, analysis, sequtils, peptutils, mhclearn
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

path = os.path.dirname(os.path.abspath(__file__))
testdir = os.path.join(path, 'testing')
datadir = os.path.join(path, 'mhcdata')

class PredictorTests(unittest.TestCase):
    """Basic tests for predictor"""

    def setUp(self):
        self.m2alleles = base.get_preset_alleles('mhc2_supertypes')
        self.peptides = peptutils.create_random_sequences(50)
        self.genbankfile = os.path.join(testdir, 'zaire-ebolavirus.gb')
        self.fastafile = os.path.join(testdir, 'zaire-ebolavirus.faa')
        self.df = sequtils.genbank_to_dataframe(self.genbankfile, cds=True)
        self.testdir = testdir
        if not os.path.exists(self.testdir):
            os.mkdir(self.testdir)
        return

    def test_peptide_utils(self):

        s = peptutils.create_random_sequences(100)
        print (s)
        seq = 'MTDDPGSGFTTVWNAVVSELNGDPKVDDGP'
        f = peptutils.get_fragments(seq=seq)
        print (f)
        return

    def test_classes(self):

        cl = base.get_predictor_classes()
        for c in cl:
            P=base.get_predictor(c)
            print (P)

    def test_cutoffs(self):

        cl = base.get_predictor_classes()
        for c in cl:
            P=base.get_predictor(c)
            P.get_allele_cutoffs()

    def test_basicmhc1(self):

        P=base.get_predictor('basicmhc1')
        print (P)
        allele='HLA-A*01:01'
        data = mhclearn.get_training_set(allele)
        peps = list(data.peptide)
        P.predict_peptides(peps[:50],alleles=allele)
        return

    def test_tepitope(self):
        """Tepitope test"""

        df = self.df
        P = base.get_predictor('tepitope')
        alleles = ["HLA-DRB1*0101", "HLA-DRB1*0305"]
        print (P)
        P.predict_proteins(df, length=11, alleles=alleles,
                          path=self.testdir)
        P.get_binders(data=P.data)
        return

    def test_netmhcpan(self):
        """netMHCpan test"""

        #requires netmHCpan is installed
        df = self.df
        P = base.get_predictor('netmhcpan')
        print (P)
        seqs = peptutils.create_random_sequences(10)
        P.predict_peptides(seqs, alleles=['HLA-A*02:02'], threads=1)
        print (len(P.data))
        return

    '''def test_netmhciipan(self):
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
        return'''

    '''def test_iedbmhc1(self):
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
            print (P.name, m)
            P.predict_proteins(df, length=8, alleles=alleles,
                                    method=m)
            b = P.get_binders(data=P.data, cutoff=5, cutoff_method='rank')
            print ('%s binders' %len(b))
        return'''

    def test_peptide_prediction(self):

        m2alleles = base.get_preset_alleles('mhc2_supertypes')
        P = base.get_predictor('tepitope')
        x = P.predict_peptides(self.peptides, alleles=self.m2alleles)
        return

    def test_multiproc(self):
        P = base.get_predictor('tepitope')
        x = P.predict_peptides(self.peptides, alleles=self.m2alleles, threads=2)
        return

    def test_fasta(self):
        """Test fasta predictions"""

        df = sequtils.fasta_to_dataframe(self.fastafile)
        alleles = ["HLA-DRB1*0101"]
        P = base.get_predictor('tepitope')
        P.predict_proteins(df, length=11, alleles=alleles, path=self.testdir)
        return

    def test_load(self):
        """Test re-loading predictions"""

        infile = os.path.join(self.testdir, 'ZEBOVgp1.csv')
        P = base.get_predictor('iedbmhc1')
        P.load(infile)
        return

    def test_features(self):
        """Test genbank feature handling"""

        df = sequtils.fasta_to_dataframe(self.fastafile)
        name = 'ZEBOVgp1'
        sequtils.dataframe_to_fasta(df)
        sequtils.check_tags(df)
        return

    def test_mhcflurry(self):
        """Test mhcflurry predictor"""

        P = base.get_predictor('mhcflurry')
        print (P)
        seqs = peptutils.create_random_sequences(10)
        P.predict_peptides(seqs, alleles=['HLA-A*02:02'], threads=1)
        print (len(P.data))
        return

    def test_iedbmhc1(self):
        """iedbmhc1 test"""

        df = self.df
        P = base.get_predictor('iedbmhc1')
        if P.check_install() == 0:
            return
        seqs = peptutils.create_random_sequences(10, length=11)
        P.predict_peptides(seqs, alleles=['HLA-A*02:02'], threads=1)
        print (len(P.data))
        return

    def quit(self):
        self.app.quit()

def run():
    unittest.main()

if __name__ == '__main__':
    unittest.main()
