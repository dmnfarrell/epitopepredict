#!/usr/bin/env python

"""
    Module that implements the TEPITOPEPan method. Includes methods for pickpocket
    and pseudosequence similarity calcaulation.
    References:
    [1] L. Zhang, Y. Chen, H.-S. Wong, S. Zhou, H. Mamitsuka, and S. Zhu, "TEPITOPEpan: extending
    TEPITOPE for peptide binding prediction covering over 700 HLA-DR molecules.,"
    PLoS One, vol. 7, no. 2, p. e30483, Jan. 2012.
    [2] H. Zhang, O. Lund, and M. Nielsen, "The PickPocket method for predicting binding
    specificities for receptors based on receptor pocket similarities: application to MHC-peptide binding."
    Bioinformatics, vol. 25, no. 10, pp. 1293-9, May 2009.
    Created January 2014
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import os, string, csv, glob
import time
from operator import itemgetter
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import PDB
from . import peptutils, utilities

refalleles = ['HLA-DRB1*0101','HLA-DRB1*0301','HLA-DRB1*0401','HLA-DRB1*0402',
           'HLA-DRB1*0404', 'HLA-DRB1*0701','HLA-DRB1*0801','HLA-DRB1*1101',
           'HLA-DRB1*1302', 'HLA-DRB1*1501','HLA-DRB5*0101']
bola = ['BoLA-DRB3*2005','BoLA-DRB3*1601','BoLA-DRB3*3301',
        'BoLA-DRB3*4801','BoLA-DRB3*4701','BoLA-DRB3*6701','BoLA-DRB3*3701']

path = os.path.dirname(os.path.abspath(__file__)) #path to module
tepitopedir = os.path.join(path,'tepitope')
home = os.path.expanduser("~")
datadir = os.path.join(path, 'mhcdata')
expdatadir = os.path.join(datadir, 'expdata')
aa_codes3 = {'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
        'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR', \
        'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA', \
        'G':'GLY', 'P':'PRO', 'C':'CYS'}

sim_matrices = ['blosum50','blosum62','pmbec']
alpha = 10
#pseudo_residues = sorted(set([j for i in pp.values() for j in i]))
pseudo_residues = [9,11,13,26,28,30,37,47,57,60,61,67,70,71,74,77,78,81,82,85,86,89]
drb_aln_file = os.path.join(tepitopedir,'bola_hla.drb.txt')

def get_matrix(name):
    if name not in sim_matrices:
        print ('no such matrix')
        return
    m = pd.read_csv(os.path.join(datadir, '%s.csv' %name),index_col=0)
    return m

blosum62 = get_matrix('blosum62')

def get_pocket_positions():
    cr = csv.reader(open(os.path.join(tepitopedir, 'tepitope_pockets.txt')))
    d = {}
    for r in cr:
        d[int(r[0])] = [int(i) for i in r[1:]]
    return d

def generate_pssm(expdata):
    """Create pssm for known binding data given a set
       of n-mers and binding score"""
    return

def get_pssms():
    """Get tepitope pssm data"""

    path = os.path.join(tepitopedir, 'pssm')
    pssms = {}
    for f in glob.glob(os.path.join(path,'*.csv')):
        name = os.path.splitext(os.path.basename(f))[0]
        pssm = pd.read_csv(f, index_col=0)
        pssm.columns = range(1,10)
        #print pssm
        pssms[name] = pssm
    return pssms

def get_pssm_score(seq, pssm):
    """Get sequence score for a given pssm"""

    total=0
    for pos in range(len(seq)):
        aa = seq[pos]
        ind = pos+1
        if aa not in pssm:
            continue
        val = pssm[aa][ind]
        if val != '-' and val != 0:
            total += val
    if total < -10:
        total = -10
    return total

def score_peptide(seq, pssm):
    """Score a single sequence in 9-mer frames"""

    nmers, s = peptutils.create_fragments(seq=seq, length=9)
    scores=[]
    for f in nmers:
        sc = get_pssm_score(f, pssm)
        pos = nmers.index(f)
        scores.append((f,pos,sc))
        #print f, sc
    return scores

def get_scores(pssm, sequence=None, peptides=None, length=11, overlap=1):
    """Score multiple fragments of a sequence in seperate fragments"""

    if peptides is None:
        peptides, s = peptutils.create_fragments(seq=sequence, length=length, overlap=overlap)
    scores=[]
    pos=0
    for p in peptides:
        sc = score_peptide(p, pssm) #get best 9-mer core
        sc = sorted(sc, key=itemgetter(2),reverse=True)
        #print (sc)
        core,i,best = sc[0]
        #print (p,core,pos,best)
        scores.append((p,core,pos,best))
        pos+=overlap
    return scores

def get_pseudo_sequence(query, positions=None, offset=28):
    """Get non redundant pseudo-sequence for a query. Assumes input is a
       sequence from alignment of MHC genes.
       """

    seq = []
    if positions == None:
        positions = pseudo_residues

    for i in positions:
        idx = i+offset
        if idx>=len(query):
            s='-'
        else:
            s=query[idx]
        #print i,s
        seq.append(s)
    return seq

def get_pockets_pseudo_sequence(query, offset=28):
    """Get pockets pseudo-seq from sequence and pocket residues.
    Args:
        query: query sequence
        offset: seq numbering offset of alignment numbering to pickpocket
        residue values
    """

    pp = pocket_residues
    seq = []
    for pos in pp:
        s=''
        for i in pp[pos]:
            s+=query[i+offset]
        seq.append(s)
    return seq

def get_allele_pocket_sequences(allele):
    """Convenience for getting an allele pocket aas"""

    alnindex = dict([(a.id,a) for a in drbaln])
    ref = alnindex[allele]
    return get_pockets_pseudo_sequence(ref)

def convert_allele_names(seqfile):
    """Convert long IPD names to common form.
    Args:
        fasta sequence file
    Returns:
        new list of seqrecords
    """

    recs = list(SeqIO.parse(seqfile,'fasta'))
    new = []
    found=[]
    for r in recs:
        a = r.description.split()[1][:10]
        a = 'HLA-'+a.replace(':','')
        if not a in found:
            found.append(a)
            s = SeqRecord(r.seq, id=a, description=r.description)
            new.append(s)
            #print (a, r.description)
    print ('%s sequences converted' %len(recs))
    filename = 'convertednames.fa'
    SeqIO.write(new, filename, 'fasta')
    return recs

def similarity_score(matrix, ref, query):
    """
    Similarity for pseudosequences using a substitution matrix.
    Args:
        matrix: subs matrix as dictionary
        ref: reference sequence
        query: query sequence
    Returns:
        a similarity value normalized to matrix
    """

    r=ref; q=query
    sim = sum([matrix[i][j] for i,j in zip(r,q) if (i!= '-' and j!='-')])
    sim1 = sum([matrix[i][j] for i,j in zip(r,r) if (i!= '-' and j!='-')])
    sim2 = sum([matrix[i][j] for i,j in zip(q,q) if (i!= '-' and j!='-')])
    #normalise the score
    normsim = sim / np.sqrt(sim1 * sim2)
    return normsim

def pickpocket(pos, allele):
    """Derive weights for a query allele using pickpocket method. This uses the
     pocket pseudosequences to determine similarity to the reference. This relies on
     the DRB alignment present in the tepitope folder.
    Args:
        pos: pocket position
        allele: query allele
    Returns:
        set of weights for library alleles at this position
    """

    alnindex = dict([(a.id,a) for a in drbaln])
    if allele not in alnindex:
        #print ('no such allele')
        return
    pos=pos-1
    #get query pseudosequence
    query = alnindex[allele]
    qp = get_pockets_pseudo_sequence(query)[pos]
    #derive weight per libary allele using similarity measure
    S = {}
    for k in librarypssms.keys():
        ref = alnindex[k]
        rp = get_pockets_pseudo_sequence(ref)[pos]
        sim = similarity_score(blosum62, rp, qp)
        S[k] = sim
        #print ind, qp, rp, ref.id, sim

    total = sum([np.power(S[k],alpha) for k in S])
    weights = dict([(k,round(np.power(S[k],alpha)/total,3)) for k in S])
    #print weights
    return weights

def create_virtual_pssm(allele):
    """Create virtual matrix from pickpocket profile weights"""

    lpssms = librarypssms
    ignore = [5,8]
    M=[]
    for i in pocket_residues:
        w = pickpocket(i, allele)
        if w == None: return
        if i in ignore:
            v = pd.Series([lpssms[l][i] for l in lpssms][0])
            M.append(v)
            continue
        #get weighted average over ref alleles for this pocket
        #and put in a dataframe
        v = pd.DataFrame([lpssms[l][i] * w[l] for l in lpssms])
        vq = v.sum()
        vq.name=i
        vq[vq<-100]=-999
        #sum all the weighted contributions
        M.append(vq)

    result = pd.DataFrame(M)
    result = result.transpose()
    #print result
    #result.to_csv(allele+'_pssm.csv',float_format='%.3f')
    return result

def allelenumber(x):
    return int(x.split('*')[1])

def get_alleles():
    """Get all alleles covered by this method."""

    df=pd.read_csv(os.path.join(tepitopedir ,'alleles.txt'))
    a= df['allele']
    return list(a)

def _get_bola_alleles():
    ref='HLA-DRB1*0101'
    ids = [2005, 1601, 3301, 4801, 4701, 6701, 3701, 2101,
           2002, 3001, 1901, 1902, 4901, 4101, 4802, 2003, 1602, 3801]
    alleles = ['BoLA-DRB3*%s' %i for i in ids]
    alleles.append(ref)
    return alleles

def get_similarities(allele, refalleles, alnindex, matrix):
    """Get distances between a query and set of ref pseudo-seqs"""

    query = alnindex[allele]
    qp = ''.join(get_pseudo_sequence(query))
    sims = []
    #for k in librarypssms.keys():
    for k in refalleles:
        ref = alnindex[k]
        rp = ''.join(get_pseudo_sequence(ref))
        #print (qp,rp)
        sim = similarity_score(matrix, rp, qp)
        sims.append((k,sim))
    return sims, qp

def compare_tepitope_alleles(alnindex):
    """Compare a set of alleles to Tepitope library HLAs"""

    t = pd.read_csv(os.path.join(tepitopedir, 'tepitope_alleles.csv'))
    alleles = t.name[:10]
    refalleles = librarypssms.keys()
    df = compare_alleles(alleles, refalleles, alnindex, reduced=False)
    return df

def compare_alleles(alleles1, alleles2, alnindex, reduced=True, cutoff=.25,
                    matrix=None, matrix_name='blosum62'):
    """Compare 2 sets of alleles for pseudo-seq distances"""

    matrix = get_matrix(matrix_name)
    #print (matrix)
    data=[]
    pseqs = {}
    if reduced==True:
        alleles1 = reduce_alleles(alleles1)
        alleles2 = reduce_alleles(alleles2)
    for a in alleles2:
        d,qp = get_similarities(a,alleles1,alnindex, matrix=matrix)
        d = pd.DataFrame(d,columns=['ref',a])
        #print (d)
        d.set_index('ref',inplace=True)
        data.append(d)
        pseqs[a]=qp

    df = pd.concat(data,axis=1)
    df = df.apply(lambda x: 1-x)
    df = df.transpose()
    df = df.sort_index()
    df['mean'] = df.mean(axis=1).round(2)
    df['nearest'] = df.min(axis=1).round(2)
    df.sort_values(['mean'], inplace=True)
    bins = np.linspace(0, 0.7, 30)

    print
    #print ('most similar alleles:')
    h = df[df['nearest']<cutoff]
    #print (h)
    h = h.drop(['mean','nearest'],axis=1)
    h = h.reindex(h.mean().sort_values().index, axis=1)

    return h

def compare(file1, file2, alnindex, reduced=True):
    """All vs all for 2 sets of sequence files"""

    hlagrps = pd.read_csv('hla-dr_groups.csv')
    #alleles1 = hlagrps.allele
    recs1 = list(SeqIO.parse(file1,'fasta'))
    recs2 = list(SeqIO.parse(file2,'fasta'))
    alleles1 = [rec.id for rec in recs1]
    alleles2 = [rec.id for rec in recs2]
    df = compare_alleles(alleles1, alleles2, alnindex, reduced)
    return df

def reduce_alleles(alleles):
    """Reduce alleles to repr set based on names"""
    red={}
    for a in alleles:
        r = a[:-2]
        if r not in red: red[r] = a
    #print red
    return red.values()

def benchmark():
    pssms = get_pssms()
    m = pssms['HLA-DRB1*0101']
    exp = utilities.getKnownMHCIIStructures()
    for i in exp:
        nativecore = exp[i]['Core']
        peptide = exp[i]['Peptide']
        allele = exp[i]['Allele']
        if allele in pssms:
            sc = get_scores(m, peptide)
            m = pssms[allele]
            m = m.transpose().to_dict()
            sc = sorted(sc, key=itemgetter(1),reverse=True)
            print (i,allele,peptide,nativecore,sc[0])

    allele = 'HLA-DRB1*0101'
    m = pssms[allele]
    m = m.transpose().to_dict()
    expfile = 'peptide_affinity_dataset/%s.txt' %allele
    exp = pd.read_csv(os.path.join(expdatadir,expfile))
    seqs = exp['SEQ_SEQUENCE']
    aff = exp['COMP_IC50']
    for s,val in zip(seqs,aff):
        sc = get_scores(m, s)
        sc = sorted(sc, key=itemgetter(1),reverse=True)
        #print s, sc[0],val
    return

def show_pocket_residues(pdbfile):
    """Test to show the pocket residues in a pdb structure"""

    from . import pymolutils
    res = pseudo_residues
    pymolutils.show_protein(pdbfile, save=True)
    pymolutils.show_residues(coords=res,offset=0)
    pymolutils.save('pockets.pse')
    return

def compare_ref(query1, query2, ref, alnindex):
    """Compare different alleles distances to reference"""

    df1=compare(ref, query1, alnindex, reduced=False)
    df2=compare(ref, query2, alnindex, reduced=False)
    bins = np.linspace(0, 0.7, 30)
    f=plt.figure()
    ax=f.add_subplot(111)
    ax.hist(df1.nearest, bins, alpha=0.5, normed=True)
    ax.hist(df2.nearest, bins, alpha=0.5, normed=True)
    ax.legend(['TEPITOPEpan HLA-DR alleles','BoLA-DRB3 alleles'])
    ax.set_xlabel('nearest neighbour distance')
    ax.set_ylabel('normalised count')
    plt.tight_layout()
    return

def test():

    aln = drbaln
    alnindex = dict([(a.id,a) for a in aln])
    compare_tepitope_alleles(alnindex)
    #d1 = compare(ref1, ref2, alnindex)
    #x = d1.merge(d2,right_index=1,left_index=1)
    #print len(x)
    #compare_ref(hla,bola,ref,alnindex)
    plt.show()
    return

pocket_residues = get_pocket_positions()
librarypssms = get_pssms()
#drb MHC alignment using IPD sequences, includes BoLA-DRB3 sequences
drbaln = AlignIO.read(drb_aln_file, "fasta")

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-t", "--test", dest="test", action='store_true',
                            help="test")
    opts, remainder = parser.parse_args()
    if opts.test == True:
        test()

if __name__ == '__main__':
    main()
