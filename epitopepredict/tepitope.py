#!/usr/bin/env python

"""
    Module that implements the TEPITOPEPan method
    Created January 2014
    Copyright (C) Damien Farrell
"""

import os, string, csv, glob
import time
from operator import itemgetter
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import PDB
import peptides, utilities

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
AAcodes3 = {'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
        'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR', \
        'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA', \
        'G':'GLY', 'P':'PRO', 'C':'CYS'}
blosum62 = pd.read_csv(os.path.join(datadir, 'blosum62.csv'),index_col=0)
blosum50 = pd.read_csv(os.path.join(datadir, 'blosum50.csv'),index_col=0)
alpha = 10

def getPocketPositions():
    cr = csv.reader(open(os.path.join(tepitopedir, 'tepitope_pockets.txt')))
    d = {}
    for r in cr:
        d[int(r[0])] = [int(i) for i in r[1:]]
    return d

def generatePSSM(expdata):
    """Create pssm for known binding data given a set
       of n-mers and binding score"""

    return

def getPSSMs():
    path = os.path.join(tepitopedir, 'pssm')
    pssms = {}
    for f in glob.glob(os.path.join(path,'*.csv')):
        name = os.path.splitext(os.path.basename(f))[0]
        pssm = pd.read_csv(f, index_col=0)
        pssm.columns = range(1,10)
        #print pssm
        pssms[name] = pssm
    return pssms

def getPSSMScore(seq, pssm):
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

def scorePeptide(seq, pssm):
    """Score a single sequence in 9-mer frames"""

    nmers, s = peptides.createFragments(seq=seq, length=9)
    scores=[]
    for f in nmers:
        sc = getPSSMScore(f, pssm)
        pos = nmers.index(f)
        scores.append((f,pos,sc))
        #print f, sc
    return scores

def getScores(pssm, sequence=None, peptide=None, length=11):
    """Score multiple fragments of a sequence in seperate fragments"""

    if peptide == None:
        peptide, s = peptides.createFragments(seq=sequence, length=length)
    scores=[]
    pos=0
    for p in peptide:
        sc = scorePeptide(p, pssm) #get best 9-mer core
        sc = sorted(sc, key=itemgetter(2),reverse=True)
        core,i,best = sc[0]
        #print (p,core,pos,best)
        scores.append((p,core,pos,best))
        pos+=1
    return scores

def getPseudoSequence(pp, query, method='tepitope'):
    """Get non redundant pseudo-seq"""

    offset=28
    seq = []
    if method == 'tepitope':
        res = pseudores
    elif method == 'netmhciipan':
        res = [9,11,13,14,26,28,30,47,56,57,60,67,70,71,74,77,78,81,85,86,89]
    for i in res:
        s=''
        s+=query[i+offset]
        #print i,s
        seq.append(s)
    return seq

def getPocketsPseudoSequence(pp, query):
    """Get pockets pseudo-seq from alignment and pocket residues"""

    offset=28 #seq numbering offset in IPD_MHC aligments
    seq = []
    for pos in pp:
        s=''
        for i in pp[pos]:
            s+=query[i+offset]
        seq.append(s)
    return seq

def getAllelePocketSequences(allele):
    """Convenience for getting an allele pocket aas"""
    alnindex = dict([(a.id,a) for a in drbaln])
    ref = alnindex[allele]
    return getPocketsPseudoSequence(pp,ref)

def convertAlleleNames(seqfile):
    """Convert long IPD names to common form"""

    recs = list(SeqIO.parse(seqfile,'fasta'))
    new = []
    found=[]
    for r in recs:
        a = r.description.split()[1][:10]
        a = 'HLA-'+a.replace(':','')
        if not a in found:
            found.append(a)
            s = SeqRecord(r.seq, id=a, description='')
            new.append(s)
            print a, r.description
    filename = 'convertednames.fa'
    SeqIO.write(new, filename, 'fasta')
    return filename

def getAlignments(file1, file2):
    """Align multiple sequences from 2 fasta files"""

    recs1 = list(SeqIO.parse(file1,'fasta'))
    recs2 = list(SeqIO.parse(file2,'fasta'))
    allrecs = recs2 + recs1
    alnfile = 'queryaln.fa'
    SeqIO.write(allrecs, alnfile, 'fasta')
    print 'doing multiple sequence alignment for %s recs' %len(allrecs)
    aln = Genome.muscleAlignment(alnfile)
    return aln

def similarityScore(blosum, r, q):
    """Similarity for pseudosequences using blosum matrix"""

    sim = sum([blosum[i][j] for i,j in zip(r,q) if (i!= '-' and j!='-')])
    sim1 = sum([blosum[i][j] for i,j in zip(r,r) if (i!= '-' and j!='-')])
    sim2 = sum([blosum[i][j] for i,j in zip(q,q) if (i!= '-' and j!='-')])
    normsim = sim / np.sqrt(sim1 * sim2)
    return normsim

def pickpocket(ind, allele):
    """Derive weights for a query allele using pickpocket method"""

    alnindex = dict([(a.id,a) for a in drbaln])
    if allele not in alnindex:
        #print 'no such allele'
        return
    ind=ind-1
    #get query pseudosequence
    query = alnindex[allele]
    qp = getPocketsPseudoSequence(pp,query)[ind]
    #derive weight per libary allele using similarity measure
    S = {}
    for k in librarypssms.keys():
        ref = alnindex[k]
        rp = getPocketsPseudoSequence(pp, ref)[ind]
        sim = similarityScore(blosum62, rp, qp)
        S[k] = sim
        #print ind, qp, rp, ref.id, sim

    total = sum([np.power(S[k],alpha) for k in S])
    weights = dict([(k,round(np.power(S[k],alpha)/total,3)) for k in S])
    #print weights
    return weights

def createVirtualPSSM(allele):
    """Create virtual matrix from pickpocket profile weights"""

    lpssms = librarypssms
    ignore = [5,8]
    M=[]
    for i in pp:
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

def getAlleles():
    """Get all alleles covered """
    df=pd.read_csv(os.path.join(tepitopedir ,'alleles.txt'))
    a= df['allele']
    return list(a)

def getBolaAlleles():
    ref='HLA-DRB1*0101'
    ids = [2005, 1601, 3301, 4801, 4701, 6701, 3701, 2101,
           2002, 3001, 1901, 1902, 4901, 4101, 4802, 2003, 1602, 3801]
    alleles = ['BoLA-DRB3*%s' %i for i in ids]
    alleles.append(ref)
    return alleles

def getSimilarities(allele, refalleles, alnindex, matrix=blosum62):
    """Get distances between a query and set of ref pseudo-seqs"""

    query = alnindex[allele]
    #qp = ''.join(getPocketsPseudoSequence(pp,query))
    qp = ''.join(getPseudoSequence(pp,query))
    sims = []
    #for k in librarypssms.keys():
    for k in refalleles:
        ref = alnindex[k]
        #rp = ''.join(getPocketsPseudoSequence(pp, ref))
        rp = ''.join(getPseudoSequence(pp, ref))
        #print qp,rp
        sim = similarityScore(matrix, rp, qp)
        sims.append((k,sim))
    return sims, qp

def plothists(df1, df2, bins=25):

    f=plt.figure()
    ax=f.add_subplot(111)
    ax.hist(df1.nearest, bins, alpha=0.5, normed=True)
    ax.hist(df2.nearest, bins, alpha=0.5, normed=True)
    plt.tight_layout()
    return ax

def plotheatmap(df):

    fig=plt.figure()
    ax=fig.add_subplot(111)
    hm=ax.pcolor(df,cmap='RdBu')
    fig.colorbar(hm, ax=ax)
    ax.set_xticks(np.arange(0.5, len(df.columns)))
    ax.set_yticks(np.arange(0.5, len(df.index)))
    ax.set_xticklabels(df.columns, minor=False, fontsize=10,rotation=75)
    ax.set_yticklabels(df.index, minor=False,fontsize=10)
    hm.set_clim(0,1)
    #from matplotlib.patches import Rectangle
    #ax.add_patch(Rectangle((0, 0), 9,ax.set_ylim()[1]-1,fill=0,lw=2,edgecolor='black',alpha=0.8))
    plt.tight_layout()
    return

def compareTepitopeAlleles(alnindex):
    """Compare a set of alleles to Tepitope library HLAs"""

    t = pd.read_csv(os.path.join(tepitopedir, 'tepitope_alleles.csv'))
    alleles = t.name[:10]
    refalleles = librarypssms.keys()
    df = compareAlleles(alleles, refalleles, alnindex, reduced=False)
    return df

def compareAlleles(alleles1, alleles2, alnindex, reduced=True):
    """Compare 2 sets of alleles for pseudo-seq distances"""

    data=[]
    pseqs = {}
    if reduced==True:
        alleles1 = reduceAlleles(alleles1)
        alleles2 = reduceAlleles(alleles2)
    for a in alleles2:
        d,qp = getSimilarities(a,alleles1,alnindex)
        d = pd.DataFrame(d,columns=['ref',a])
        d.set_index('ref',inplace=True)
        data.append(d)
        pseqs[a]=qp

    df = pd.concat(data,axis=2)
    df = df.apply(lambda x: 1-x)
    df = df.transpose()
    df = df.sort_index()
    df['mean'] = df.mean(axis=1).round(2)
    df['nearest'] = df.min(axis=1).round(2)
    df.sort(['nearest'], inplace=True)
    bins=np.linspace(0, 0.7, 30)
    df.hist(column=['nearest'],bins=bins,grid=0,color='gray')
    df.to_csv('allele_similarities.csv')
    #plt.suptitle('bola-drb3 pseudo-sequence distances')
    #plt.savefig('allele_sims_hist.png')
    #plt.show()
    #plt.clf()
    print
    print 'most similar alleles:'
    h = df[df['nearest']<0.25]
    print h[['nearest','mean']].sort()
    h = h.drop(['mean','nearest'],axis=1)
    h = h.reindex_axis(h.mean().order().index, axis=1)
    plotheatmap(h)
    found = list(df.index)
    #print found
    for r in refalleles:
        pseqs[r] =  ''.join(getPseudoSequence(pp, alnindex[r]))
        if r not in found:
            found.append(r)
    for i in sorted(pseqs):
        print '%-15s' %i, pseqs[i]
    #distanceTree(seqs=[SeqRecord(Seq(pseqs[i]),i) for i in found], ref=refalleles[0])
    #ETETree(seqs=[SeqRecord(Seq(pseqs[i]),i) for i in found],
    #        ref=refalleles[0],metric=dict(df['nearest']))
    return h

def compare(file1, file2, alnindex, reduced=True):
    """All vs all for 2 sets of sequence files"""

    hlagrps = pd.read_csv('hla-dr_groups.csv')
    #alleles1 = hlagrps.allele
    recs1 = list(SeqIO.parse(file1,'fasta'))
    recs2 = list(SeqIO.parse(file2,'fasta'))
    alleles1 = [rec.id for rec in recs1]
    alleles2 = [rec.id for rec in recs2]
    df = compareAlleles(alleles1, alleles2, alnindex, reduced)
    return df

def reduceAlleles(alleles):
    """Reduce alleles to repr set based on names"""
    red={}
    for a in alleles:
        r = a[:-2]
        if r not in red: red[r] = a
    #print red
    return red.values()

def benchmark():
    pssms = getPSSMs()
    m = pssms['HLA-DRB1*0101']
    exp = utilities.getKnownMHCIIStructures()
    for i in exp:
        nativecore = exp[i]['Core']
        peptide = exp[i]['Peptide']
        allele = exp[i]['Allele']
        if allele in pssms:
            sc = getScores(m, peptide)
            m = pssms[allele]
            m = m.transpose().to_dict()
            sc = sorted(sc, key=itemgetter(1),reverse=True)
            print i,allele,peptide,nativecore,sc[0]

    allele = 'HLA-DRB1*0101'
    m = pssms[allele]
    m = m.transpose().to_dict()
    expfile = 'peptide_affinity_dataset/%s.txt' %allele
    exp = pd.read_csv(os.path.join(expdatadir,expfile))
    seqs = exp['SEQ_SEQUENCE']
    aff = exp['COMP_IC50']
    for s,val in zip(seqs,aff):
        sc = getScores(m, s)
        sc = sorted(sc, key=itemgetter(1),reverse=True)
        #print s, sc[0],val
    return

def compareRef(query1, query2, ref, alnindex):
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

    #alnindex = dict([(a.id,a) for a in aln])
    #Genome.showAlignment(aln)
    aln = drbaln
    alnindex = dict([(a.id,a) for a in aln])
    compareTepitopeAlleles(alnindex)
    #d1 = compare(ref1, ref2, alnindex)
    #x = d1.merge(d2,right_index=1,left_index=1)
    #print len(x)
    #compareRef(hla,bola,ref,alnindex)
    plt.show()
    return

#declare these as global for convenience
pp = getPocketPositions()
pseudores = sorted(set([j for i in pp.values() for j in i]))
librarypssms = getPSSMs()
#drb HLA + BOLA alignments
drbaln = AlignIO.read(os.path.join(tepitopedir,'bola_hla.drb.txt'), "fasta")

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
