#!/usr/bin/env python

"""
    Module implementing peptide sequence/structure utilities.
    Created March 2013
    Copyright (C) Damien Farrell
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

from __future__ import absolute_import, print_function
import os, random, csv
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
#from Bio.Alphabet import IUPAC
from . import utilities

AAletters = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P',\
             'S', 'R', 'T', 'W', 'V', 'Y']
AAcodes3 ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR', \
'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA', \
'G':'GLY', 'P':'PRO', 'C':'CYS'}

def create_random_sequences(size=100,length=9):
    """Create library of all possible peptides given length"""

    aa = AAletters
    vals=[]
    for i in range(0, size):
        seq = ""
        for s in range(0,length):
            index = random.randint(0, len(aa)-1)
            seq += aa[index]
        vals.append(seq)
    return vals

def create_random_peptides(size=100,length=9):
    """Create random peptide structures of given length"""

    seqs = create_random_sequences(size, length)
    files=[]
    for seq in seqs:
        sfile = buildPeptideStructure(seq)
        files.append(sfile)
    return files

def get_fragments(seq=None, overlap=1, length=11, **kwargs):
    """Generate peptide fragments from a sequence.
    Returns:
        dataframe of peptides with position column.
    """

    frags=[]
    for i in range(0,len(seq),overlap):
        if i+length>len(seq): continue
        frags.append([i+1,seq[i:i+length]])
    df = pd.DataFrame(frags, columns=['pos','peptide'])
    return df

def create_fragments(protfile=None, seq=None, length=9, overlap=1, quiet=True):
    """generate peptide fragments from a sequence"""

    if protfile != None:
        rec = SeqIO.parse(open(protfile,'r'),'fasta').next()
        seq = str(rec.seq)
    frags=[]
    for i in range(0,len(seq),overlap):
        if i+length>len(seq): continue
        frags.append(seq[i:i+length])
    if quiet == False:
        print ('created %s fragments of length %s' %(len(frags),length))
    if protfile != None:
        cw =  csv.writer(open('%s_fragments' %protfile,'w'))
        for f in frags:
            cw.writerow([f])
    return frags, seq

def get_all_fragments(exp, length=11):

    P=pd.read_csv(exp)
    peptides = P['SEQ_SEQUENCE']
    allseqs=[]
    for p in peptides:
        seqs, seqstr = createFragments(seq=p,length=length)
        allseqs.extend(seqs)
    cw=csv.writer(open('all_fragments.csv','w'))
    for i in allseqs:
        cw.writerow([i])
    print ('got %s fragments' %len(allseqs))
    return allseqs

def get_AAsubstitutions(template):
    """
    Get all the possible sequences from substituting every AA
      into the given sequence at each position. This gives a total of
       19 by n amino acid positions.
    """

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

def get_AAfraction(seq, amino_acids=None):
    """Get fraction of give amino acids in a sequence"""

    X = ProteinAnalysis(seq)
    #h = X.protein_scale(ProtParam.ProtParamData.kd, len(seq), 0.4)
    nonpolar = ['A','V','L','F','I','W','P']
    if amino_acids == None:
        amino_acids = nonpolar
    count=0
    for aa, i in X.count_amino_acids().items():
        if aa in amino_acids:
            count+=i
    if count == 0: return 0
    frac = round(float(count)/len(seq),2)
    return frac

def net_charge(seq):
    """Get net charge of a peptide sequence"""

    X = ProteinAnalysis(seq)
    ac = 0
    ba = 0
    for aa, i in X.count_amino_acids().items():
        if aa in ['D','E']:
            ac -= i
        elif aa in ['K','R']:
            ba += i
    return ac + ba

def compare_anchor_positions(x1, x2):
    """Check if anchor positions in 9-mers are mutated"""

    if x1 is None or x2 is None:
        return 0
    p1 = list(get_fragments(x1, length=9).peptide)
    p2 = list(get_fragments(x2, length=9).peptide)
    #is mutation in anchor residue
    for a,b in zip(p1,p2):
        if a[1] != b[1] or a[8] != b[8]:
            return 1
    return 0

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
