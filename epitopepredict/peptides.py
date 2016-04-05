#!/usr/bin/env python

"""
    Module implementing peptide structure utilities.
    Created March 2013
    Copyright (C) Damien Farrell
"""

import os, random, csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Alphabet import IUPAC
import utilities

AAletters = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P',\
             'S', 'R', 'T', 'W', 'V', 'Y']
AAcodes3 ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR', \
'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA', \
'G':'GLY', 'P':'PRO', 'C':'CYS'}

def createRandomSequences(size=100,length=9):
    """Create library of all possible peptides given length"""
    aa = list(IUPAC.IUPACProtein.letters)
    vals=[]
    for i in range(0, size):
        seq = ""
        for s in range(0,length):
            index = random.randint(0, len(aa)-1)
            seq += aa[index]
        vals.append(seq)
    return vals

def createRandomPeptides(size=100,length=9):
    """Create random peptide structures of given length"""
    seqs = createRandomSequences(size, length)
    files=[]
    for seq in seqs:
        sfile = buildPeptideStructure(seq)
        files.append(sfile)
    return files

def createFragments(protfile=None, seq=None, length=9, overlap=1, quiet=True):
    """generate peptide fragments from a sequence"""
    if protfile != None:
        rec = SeqIO.parse(open(protfile,'r'),'fasta').next()
        seq = str(rec.seq)
    frags=[]
    for i in range(0,len(seq),overlap):
        if i+length>len(seq): continue
        frags.append(seq[i:i+length])
    if quiet == False:
        print 'created %s fragments of length %s' %(len(frags),length)
    if protfile != None:
        cw =  csv.writer(open('%s_fragments' %protfile,'w'))
        for f in frags:
            cw.writerow([f])
    return frags, seq

def createProteinFragments(protfile, seqrange='', length=9):
    """Create peptide fragments for a protein"""
    #homologyModel()
    seqs, seqstr = createFragments(protfile, length=length)
    filenames = []
    if seqrange != '':
        lims = utilities.getListfromConfig(seqrange)
        seqrange = seqs[lims[0]:lims[1]]
    else:
        seqrange = seqs
    for s in seqrange:
        fname = buildPeptideStructure(s)
        filenames.append(fname)
    return filenames

def getAllFragments(exp, length=11):
    import pandas as pd
    P=pd.read_csv(exp)
    peptides = P['SEQ_SEQUENCE']
    allseqs=[]
    for p in peptides:
        seqs, seqstr = createFragments(seq=p,length=length)
        allseqs.extend(seqs)
    cw=csv.writer(open('all_fragments.csv','w'))
    for i in allseqs:
        cw.writerow([i])
    print 'got %s fragments' %len(allseqs)
    return allseqs

def getAASubstitutions(template):
    """Get all the possible sequences from substituting every AA
      into the given sequence at each position. This gives a total of
       19aa * n positions """
    aas = AAletters
    seqs = []
    matrix = []
    for pos in range(len(template)):
        s = template[pos]
        x=[]
        for a in aas:
            if a == s:
                continue
            newseq = list(template)
            newseq[pos] = a
            newseq = ''.join(newseq)
            seqs.append(newseq)
            x.append(newseq)
        matrix.append(x)
    return seqs, matrix

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-t", "--test", dest="test", action='store_true',
                            help="test")
    opts, remainder = parser.parse_args()

    if opts.test == True:
        getAllFragments(opts.file,9)

if __name__ == '__main__':
    main()
