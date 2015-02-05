#!/usr/bin/env python

"""
    Genome annotation methods
    Created November 2013
    Copyright (C) Damien Farrell
"""

import sys, os, shutil, string, types
import csv, glob, pickle, operator
import time, re
from collections import OrderedDict
import numpy as np
import pylab as plt
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import utilities

featurekeys = ['type','protein_id','locus_tag','gene','db_xref',
               'product', 'note', 'translation','pseudo','start','end']
typecolors = ['blue','green','brown','orange','purple','lightblue','yellow','red']


def drawFeatures(tracks, highlight=[], filename=None,
                color=None, pagesize=(600,200),name=''):
    """Draw gene features in multiple tracks for general comparison"""

    from reportlab.lib import colors
    from reportlab.lib.units import cm
    from Bio.Graphics import GenomeDiagram
    gd_diagram = GenomeDiagram.Diagram("test1")
    i=0
    locs = [f.location for f in tracks[0]]
    start = min([l.start for l in locs])
    end = max([l.end for l in locs])

    for features in tracks:
        track = gd_diagram.new_track(1, greytrack=True,greytrack_labels=2,
                                    name=name,scale_fontsize=10,scale_fontangle=0,
                                    scale_largetick_interval=1500,scale_color=colors.gray,
                                    scale_largeticks=1,scale_smallticks=0,
                                    greytrack_fontsize=20)
        gd_feature_set = track.new_set()
        #features.sort(key=lambda x: abs(x.location.end-x.location.start))

        for feature in features:
            if color!=None: c=color
            else:
                c=typecolors[i]
            if feature in highlight: c='red'
            if feature.type != 'CDS': continue
            gd_feature_set.add_feature(feature, color=c, label=True,
                                        sigil="ARROW", arrowhead_length=0.5,
                                        arrowshaft_height=0.2, tracklines=1,
                                        label_size=12, label_angle=45)
        i+=1
    gd_diagram.draw(format='linear', orientation="landscape", pagesize=pagesize,
                    fragments=1, start=start, end=end)
    if filename==None:
        filename = 'genediagram'
    gd_diagram.write(filename+'.png', "PNG")
    return filename+'.png'

def drawFeatures2(tks, highlight=[], filename=None, color=None, **kwargs):
    from biograpy import Panel, tracks, features
    panel=Panel(fig_width=500,fig_height=150, fig_dpi=100,padding=10,
                grid='both')
    for feats in tks:
        test_track = tracks.BaseTrack(name = 'X',
                                       height = .2)
        for feat in feats:
            if feat.location.strand == -1:
                arr = "larrow,pad=0.6"
            else:
                arr = "rarrow,pad=0.6"
            q=feat.qualifiers
            tag=q['locus_tag'][0]
            test_track.append(features.GenericSeqFeature(feat,
                              name=tag,boxstyle=arr,
                              ))
        panel.add_track(test_track)
    if filename==None:
        filename = 'genediagram'
    panel.save(filename+'.png')
    panel.close()
    return filename+'.png'

def drawGenomeMap(infile, filename=None):
    """Draw whole circular genome"""
    from Bio.Graphics import GenomeDiagram
    from Bio.SeqUtils import GC
    from reportlab.lib import colors
    genome = SeqIO.read(infile,'genbank')
    gdd = GenomeDiagram.Diagram('test')
    gdt1 = gdd.new_track(4, greytrack=1, name='CDS', scale=0)
    gdt2 = gdd.new_track(3, greytrack=1, name='tRNA', scale=0)
    gdt3 = gdd.new_track(2, greytrack=0, name='GC content', scale=0)
    gdf1 = gdt1.new_set('feature')
    gdf2 = gdt2.new_set('feature')
    gdgs = gdt3.new_set('graph')
    graphdata = [(f.location.start,GC(f.extract(genome.seq))) for f in genome.features]
    #print graphdata
    gdgs.new_graph(graphdata, 'GC content', style='line', colour=colors.black)
    for feature in genome.features:
        if feature.type == 'CDS':
            gdf1.add_feature(feature, label=False, colour=colors.green)
        elif feature.type == 'tRNA' :
            gdf2.add_feature(feature, label=True, colour=colors.red)
    gdd.draw(format='circular', orientation='landscape',
            tracklines=0, pagesize='A4', fragments=5, circular=1)
    if filename==None:
        filename = 'genediagram.png'
    gdd.write(filename, "PNG")
    return filename

def distanceTree(seqfile=None, seqs=None, ref=None):
    """Phylo tree for sequences"""

    aln = clustalAlignment(seqfile, seqs)
    from Bio import Phylo
    tree = Phylo.read('temp.dnd', 'newick')
    leafList = tree.get_terminals()
    if ref != None:
        tree.root_with_outgroup(ref)
    f=plt.figure()
    #Phylo.draw_graphviz(tree,font_size='9', prog='neato')
    Phylo.draw(tree)
    return

def ETETree(seqs, ref, metric):
    """Tree showing bola alleles covered by tepitope"""
    from ete2 import Tree,PhyloTree,TreeStyle,NodeStyle
    aln = Genome.clustalAlignment(seqs=seqs)
    t = Tree('temp.dnd')
    #t.set_outgroup(t&ref)
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = "c"
    ts.arc_start = -180
    ts.arc_span = 180
    cutoff=0.25
    def func(node):
        if node.name=='NoName' or not node.name in metric:
            return False
        if metric[node.name]<=cutoff:
            return True
    matches = filter(func, t.traverse())
    print len(matches), "nodes have distance <=%s" %cutoff
    nst1 = NodeStyle()
    nst1["bgcolor"] = "Yellow"
    for n in matches:
        n.set_style(nst1)
    nst2 = NodeStyle()
    nst2["bgcolor"] = "LightGreen"
    hlanodes = [t.get_leaves_by_name(name=r)[0] for r in refalleles]
    for n in hlanodes:
        n.set_style(nst2)
    t.show(tree_style=ts)
    return

def doLocalBlast(database, query, output=None, maxseqs=10, evalue=0.001,
                    compress=False):
    """Blast a local db"""
    start = time.time()
    if output == None:
        output = os.path.splitext(query)[0]+'.xml'
    print output
    from Bio.Blast.Applications import NcbiblastxCommandline
    cline = NcbiblastxCommandline(query=query,  cmd='blastp', max_target_seqs=maxseqs,
                                 db=database, outfmt=5, out=output,
                                 evalue=evalue, num_threads=3)
    print cline
    stdout, stderr = cline()
    end = time.time()
    print 'BLAST done. took %s seconds' %(round(end-start,2))
    if compress == True:
        utilities.compress(output, remove=True)
    return

def parseBlastRec(rec):
    """Parse blast record alignment(s)"""
    if len(rec.alignments) == 0 : print 'no alignments'
    recs=[]
    qry = rec.query.split()[0]
    #print rec.query_length
    for align in rec.alignments:
        hsp = align.hsps[0]
        subj = align.title.split('>')[0]
        if qry == subj: continue
        recs.append([subj, hsp.score, hsp.expect, hsp.identities,
                    hsp.positives, rec.query_length,hsp.sbjct])
    return recs

def getBlastResults(handle=None, filename=None, n=80):
    """Get blast results into dataframe"""
    from Bio.Blast import NCBIXML
    import gzip
    if filename!=None:
        #handle = open(filename)
        handle = gzip.open(filename, 'rb')
    blastrecs = NCBIXML.parse(handle)
    rows=[]
    for rec in blastrecs:
        r = parseBlastRec(rec)
        rows.extend(r)
    df = pd.DataFrame(rows, columns=['subj','score','expect','identity',
                            'positive','query_length','sequence'])
    df['perc_ident'] = df.identity/df.query_length*100
    return df

def testblast():
    database = 'all_genomes'
    f='test.faa'
    #doLocalBlast(database, f, compress=True)
    df=getBlastResults(filename='test.xml.gz')
    print df[:5]
    return

def getFastaProteins(infile,idindex=0):
    """Get fasta proteins into dataframe"""

    recs = SeqIO.parse(infile,'fasta')
    keys = ['name','sequence','description']
    data = [(r.name.split('|')[idindex],r.seq.tostring(),r.description) for r in recs]
    df = pd.DataFrame(data,columns=(keys))
    df.set_index(['name'],inplace=True)
    return df

def convertSequenceFormat(infile, format='embl'):
    informat = os.path.splitext(infile)[1][1:]
    print 'input format: %s' %informat
    print 'output format: %s' %format
    count = SeqIO.convert(infile, informat, 'converted', format)
    print "Converted %i records" %count
    return

def getCDS(df):
    """Get CDS with transaltions from genome dataframe"""
    cds = df[df.type=='CDS']
    cdstrans = cds[cds.translation.notnull()]
    return cdstrans

def fastaFormatfromFeature(feature):
    """Get fasta formatted sequence from a genome feature"""
    name = feature.qualifiers['locus_tag'][0]
    if not feature.qualifiers.has_key('translation'):
        return ''
    seq = feature.qualifiers['translation'][0]
    rec = SeqRecord(Seq(seq),id=name,
                description=feature.qualifiers['product'][0])
    fastafmt = rec.format("fasta")
    return fastafmt

def getORFs(sequence, treshold):

    from Bio.Data import CodonTable
    start_codon_index = 0
    end_codon_index = 0
    start_codon_found = False
    btable = CodonTable.unambiguous_dna_by_name["Bacterial"]
    start_codons = btable.start_codons
    stop_codons = btable.stop_codons
    orfs = []

    for j in range(0, 3):
        for indx in range(j, len(sequence), 3):
            current_codon = sequence[indx:indx+3]
            if current_codon in start_codons and not start_codon_found:
                start_codon_found = True
                start_codon_index = indx
            if current_codon in stop_codons and start_codon_found:
                end_codon_index = indx
                length = end_codon_index - start_codon_index
                if length >= treshold * 3:
                    print length, start_codon_index,end_codon_index
                    orfs.append(start_codon_index)
                    if length % 3 != 0:
                        print "it's going to complain"
                    #print len(sequence)-end_codon_index-3
                    snippet = Seq(sequence[start_codon_index:end_codon_index])
                    print snippet
                    try:
                        protein = Seq.translate(snippet, table=11, cds=True)
                        print "%i %s" % (length/3, protein)
                    except Exception, e:
                        print e
                        pass
                start_codon_found = False
        start_codon_index = 0
        end_codon_index = 0
        start_codon_found = False

    return len(orfs), orfs

def dataframe2Fasta(df, seqkey='translation', idkey='locus_tag',
                     productkey='product',
                     outfile='out.faa'):
    """Genbank features to fasta file"""

    seqs=[]
    for i,row in df.iterrows():
        rec = SeqRecord(Seq(row[seqkey]),id=row[idkey],
                            description=row[productkey])
        seqs.append(rec)
    SeqIO.write(seqs, outfile, "fasta")
    return outfile

def genbank2Dataframe(infile, cds=False, quiet=True):
    """Get genome records from a genbank file into a dataframe
      returns a dataframe with a row for each cds/entry"""

    genome = SeqIO.read(infile,'genbank')
    #preprocess features
    allfeat = []
    for (item, f) in enumerate(genome.features):
        x = f.__dict__
        q = f.qualifiers
        x.update(q)
        d = {}
        d['start'] = f.location.start
        d['end'] = f.location.end
        for i in featurekeys:
            if x.has_key(i):
                if type(x[i]) is types.ListType:
                    d[i] = x[i][0]
                else:
                    d[i] = x[i]
        allfeat.append(d)

    df = pd.DataFrame(allfeat,columns=featurekeys)
    df['length'] = df.translation.str.len()
    df = checkTags(df)
    if quiet == False:
        print '---- %s summary ----' %infile
        s = genbankSummary(df)
        for i in s:
            print i,':',s[i]
    if cds == True:
        df = getCDS(df)
        df['order'] = range(1,len(df)+1)
    return df

def checkTags(df):
    """Check genbank tags to make sure they are not empty"""

    def replace(x):
        if pd.isnull(x.locus_tag):
            return x.gene
        else:
            return x.locus_tag
    df['locus_tag'] = df.apply(replace,1)
    return df

def genbankSummary(df):
    """Genbank dataframe summary"""
    def hypo(val):
        val = val.lower()
        kwds=['hypothetical','conserved protein','unknown protein']
        for k in kwds:
            if k in val:
                return True
        else:
            return False
    coding = df[df.type=='CDS']
    trna = df[df.type=='tRNA']
    products = coding[coding['product'].notnull()]
    cdstrans = coding[coding.translation.notnull()]
    hypo = products[products['product'].apply(hypo)]
    pseudo = df[(df.type=='gene') & (df.pseudo.notnull())]
    notags = df[df.locus_tag.isnull()]
    s = {}
    s['total features'] = len(df)
    s['coding sequences'] = len(coding)
    s['cds with translations'] = len(cdstrans)
    s['cds with gene names'] = len(coding.gene.dropna())
    s['hypothetical'] = len(hypo)
    s['pseudogenes'] = len(pseudo)
    s['trna'] = len(trna)
    s['no locus tags'] =  len(notags)
    if len(cdstrans)>0:
        avlength = int(np.mean([len(i) for i in cdstrans.translation]))
        s['mean sequence length'] =  avlength

    return s

def blastORFs(infile, seq):
    """Blast translations for a sequence"""
    return

def blastGenomeRecords():
    """Blast hypotheticals in genome for up to date results"""
    return

def findkeyword(f):
    """Get keyword from a field"""

    f = f[:100]
    f = re.split('[ |/,.:]+',f)
    l=[]
    for i in f:
        if i.startswith('Rv'):
            s = i.strip()
            l.append(s)
    l = list(OrderedDict.fromkeys(l))
    return l

def indexGenbankFeatures(gb_record, feature_type, qualifier):
    """Index features by qualifier value for easy access"""

    answer = dict()
    for (index, feature) in enumerate(gb_record.features):
        print index, feature
        if feature.type==feature_type:
            if qualifier in feature.qualifiers:
                values = feature.qualifiers[qualifier]
                if not type(values) is types.ListType:
                    values = [values]
                for value in values:
                    if value in answer:
                        print "WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index)
                    else:
                        answer[value] = index
    return answer

def getGenesbyLocation(genome, feature, within=20):
    """Gets all featues within a given distance of a gene"""

    start = feature.location.start
    F = []
    dists = []
    for (i, feat) in enumerate(genome.features):
        if feat.type != 'CDS':
            continue
        #print feat.location.start in feature
        dist = abs(feat.location.start - start)
        if dist < within:
            F.append((feat, dist))
            #print i, start, feat.location, feat.qualifiers['locus_tag'][0], dist
    if len(F)==0:
        return None
    F = [i[0] for i in sorted(F, key=operator.itemgetter(1))]
    return F

def getTranslation(feature, genome, cds=True):
    """Check the translation of a cds feature"""

    trtable = "Bacterial"
    q = feature.qualifiers
    #trans = q1['translation'][0]
    seq = feature.extract(genome.seq)
    e=None
    try:
        protein = seq.translate(table=trtable,cds=cds,to_stop=True)
        #print 'protein seq:',protein
    except Exception as e:
        protein = ''
    return protein, e

def pairwiseAlignment(rec1,rec2):
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62
    gap_open = -10
    gap_extend = -0.5
    alns = pairwise2.align.localds(rec1, rec2, matrix, gap_open, gap_extend)
    return alns

def clustalAlignment(filename=None, seqs=None, command="clustalw"):
    """Align 2 sequences with clustal"""

    if filename == None:
        filename = 'temp.faa'
        SeqIO.write(seqs, filename, "fasta")
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import ClustalwCommandline
    cline = ClustalwCommandline(command, infile=filename)
    #print 'performing clustal alignment..'
    stdout, stderr = cline()
    align = AlignIO.read(name+'.aln', 'clustal')
    return align

def needleAlignment(seq1,seq2,outfile='needle.txt'):
    """Align 2 sequences with needle"""

    SeqIO.write(seq1, 'alpha.faa', "fasta")
    SeqIO.write(seq2, 'beta.faa', "fasta")
    from Bio.Emboss.Applications import NeedleCommandline
    cline = NeedleCommandline(asequence='alpha.faa', bsequence='beta.faa',
                              gapopen=30, gapextend=0.5, outfile=outfile)
    stdout, stderr = cline()
    align = AlignIO.read('needle.txt',"emboss")
    return align

def muscleAlignment(filename=None, seqs=None):
    """Align 2 sequences with muscle"""

    if filename == None:
        filename = 'temp.faa'
        SeqIO.write(seqs, filename, "fasta")
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input=filename, out=name+'.txt')
    stdout, stderr = cline()
    align = AlignIO.read(name+'.txt', 'fasta')
    return align

def showAlignment(aln, diff=False, offset=0):
    """Show diff alignment"""

    ref=aln[0]
    start=0; end=80
    #offset=28
    lbls = np.arange(start,end,10)-offset
    print ('%-20s' %'name'),''.join([('%-10s' %i) for i in lbls])
    print ('%20s' %ref.id), ref.seq[start:end]
    if diff == True:
        for a in aln[1:]:
            diff=''
            for i,j in zip(ref,a):
                if i != j:
                    diff+=j
                else:
                    diff+='-'
            print ('%20s' %a.id), diff[start:end]
    else:
        for a in aln[1:]:
            print ('%20s' %a.id), a.seq[start:end]
    return

def getIdentity(aln):
    """Get sequence identity of alignment for overlapping region only"""
    j=0
    i=0
    record = aln[1]
    start=None; end=None #use these to get local overlap
    for aa in record.seq:
        aa1 = aln[0].seq[j]
        if aa == '-' or aa1 == '-':
            pass
        else:
            if aa == aa1:
                if start == None:
                    start = j
                end = j+1
                i+=1
        j+=1
    overlap = end-start
    percent = round(100.0*i/overlap,1)
    return percent, overlap

def formatAlignment(aln):
    t=''
    for i in range(0,len(aln[0]),80):
        for a in aln:
            t+=('%15s' %a.id) +' '+ a.seq.tostring()[i:i+80]+'\n'
        t+='\n'
    return t

def getAlignment(f1,f2,genome1,genome2,checkcds=True):
    p1,e1 = getTranslation(f1, genome1, checkcds)
    p2,e2 = getTranslation(f2, genome2, checkcds)
    if len(p1) == 0 or len(p2) == 0:
        return None
    aln = needleAlignment(SeqRecord(p1,'a'),SeqRecord(p2,'b'))
    return aln


def getFeatureQualifier(f, qualifier):
    if f.qualifiers.has_key(qualifier):
        fq = f.qualifiers[qualifier][0].lower()
    else:
        fq = None
    return fq


def translateSixFrame(seq):
    """Translate seq in 6 frames"""
    from cogent import DNA
    from cogent.core.genetic_code import DEFAULT as standard_code
    translations = standard_code.sixframes(seq)
    stops_frame1 = standard_code.getStopIndices(seq, start=0)
    print translations
    return


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-t", "--test", dest="test", action='store_true',
                            help="test")
    opts, remainder = parser.parse_args()

    if opts.test == True:
        testblast()


if __name__ == '__main__':
    main()
