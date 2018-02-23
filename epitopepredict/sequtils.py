#!/usr/bin/env python

"""
    Sequence utilities and genome annotation methods
    Created November 2013
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, shutil, string, types
import csv, glob, pickle, operator
import time, re
from collections import OrderedDict
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from  . import utilities

featurekeys = ['type','protein_id','locus_tag','gene','db_xref',
               'product', 'note', 'translation','pseudo','start','end']
typecolors = ['blue','green','brown','orange','purple','lightblue','yellow','red']

def draw_genome_map(infile, filename=None):
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

def distance_tree(filename=None, seqs=None, ref=None):
    """Basic phylogenetic tree for an alignment"""

    from Bio import Phylo
    if seqs is not None:
        aln = clustal_alignment(None, seqs)
        filename = 'temp.dnd'
    tree = Phylo.read(filename, 'newick')
    leaf_list = tree.get_terminals()
    if ref != None:
        tree.root_with_outgroup(ref)

    #Phylo.draw_graphviz(tree,font_size='9', prog='neato')
    f = plt.figure(figsize=(8,8))
    ax=f.add_subplot(111)
    ax.set_axis_bgcolor('white')
    Phylo.draw(tree, axes=ax)
    return tree

def ete_tree(aln):
    """Tree showing alleles"""

    from ete2 import Tree,PhyloTree,TreeStyle,NodeStyle

    t = Tree('temp.dnd')
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = "c"
    ts.arc_start = -180
    ts.arc_span = 180
    cutoff=0.25
    def func(node):
        if node.name=='NoName': #or not node.name in metric:
            return False
        #if metric[node.name]<=cutoff:
        #    return True
    matches = filter(func, t.traverse())
    print (len(matches), "nodes have distance <=%s" %cutoff)
    nst1 = NodeStyle()
    nst1["bgcolor"] = "Yellow"
    for n in matches:
        n.set_style(nst1)
    nst2 = NodeStyle()
    nst2["bgcolor"] = "LightGreen"
    #hlanodes = [t.get_leaves_by_name(name=r)[0] for r in refalleles]
    #for n in hlanodes:
    #    n.set_style(nst2)
    t.show(tree_style=ts)
    return

def local_blast(database, query, output=None, maxseqs=50, evalue=0.001,
                    compress=False, cmd='blastp', cpus=2, **kwargs):
    """Blast a local database"""

    if output == None:
        output = os.path.splitext(query)[0]+'_blast.txt'
    from Bio.Blast.Applications import NcbiblastxCommandline
    outfmt = '"6 qseqid sseqid qseq sseq pident length mismatch gapopen qstart qend sstart send evalue bitscore"'
    cline = NcbiblastxCommandline(query=query, cmd=cmd, db=database,
                                 max_target_seqs=maxseqs,
                                 outfmt=outfmt, out=output,
                                 evalue=evalue, num_threads=cpus, **kwargs)
    #print (cline)
    stdout, stderr = cline()
    return

def get_blast_results(filename):
    """
    Get blast results into dataframe. Assumes column names from local_blast method.
    Returns:
        dataframe
    """

    cols = ['qseqid','sseqid','qseq','sseq','pident','length','mismatch','gapopen',
            'qstart','qend','sstart','send','evalue','bitscore']
    res = pd.read_csv(filename, names=cols, sep='\t')
    #res = res[res['pident']>=ident]
    return res

def blast_sequences(database, seqs, labels=None, **kwargs):
    """
    Blast a set of sequences to a local blast database
    Args:
        database: local blast db name
        seqs: sequences to query, list of strings or Bio.SeqRecords
    Returns:
        pandas dataframe with top blast results
    """

    res = []
    if labels == None:
        labels = seqs
    for seq, name in zip(seqs,labels):
        if type(seq) is not SeqRecord:
            rec = SeqRecord(Seq(seq),id='temp')
        else:
            rec = seq
            name = seq.id
        SeqIO.write([rec], 'tempseq.fa', "fasta")
        local_blast(database, 'tempseq.fa', **kwargs)
        df = get_blast_results(filename='tempseq_blast.txt')
        df['id'] = name
        res.append(df)
    return pd.concat(res)

def fasta_to_dataframe(infile, header_sep=None, key='locus_tag', seqkey='translation'):
    """Get fasta proteins into dataframe"""

    recs = SeqIO.parse(infile,'fasta')
    keys = [key,seqkey,'description']
    data = [(r.name,str(r.seq),str(r.description)) for r in recs]
    df = pd.DataFrame(data,columns=(keys))
    df['type'] = 'CDS'
    #fix bad names
    if header_sep not in ['',None]:
        df[key] = df[key].apply(lambda x: x.split(header_sep)[0],1)
    df[key] = df[key].str.replace('|','_')
    return df

def convert_sequence_format(infile, outformat='embl'):
    """convert sequence files using SeqIO"""

    informat = os.path.splitext(infile)[1][1:]
    if informat == 'fa':
        informat = 'fasta'
    print ('input format: %s' %informat)
    print ('output format: %s' %outformat)
    outfile = os.path.splitext(infile)[0]+'.'+outformat
    count = SeqIO.convert(infile, informat, outfile, outformat)
    print ("Converted %i records" %count)
    return

def get_cds(df):
    """Get CDS with transaltions from genbank dataframe"""

    cds = df[df.type=='CDS']
    cdstrans = cds[cds.translation.notnull()]
    return cdstrans

def fasta_format_from_feature(feature):
    """Get fasta formatted sequence from a genome feature"""

    name = feature.qualifiers['locus_tag'][0]
    if not feature.qualifiers.has_key('translation'):
        return ''
    seq = feature.qualifiers['translation'][0]
    rec = SeqRecord(Seq(seq),id=name,
                description=feature.qualifiers['product'][0])
    fastafmt = rec.format("fasta")
    return fastafmt

def dataframe_to_seqrecords(df, seqkey='sequence', idkey='id'):
    """dataframe to list of Bio.SeqRecord objects"""

    seqs=[]
    for i,r in df.iterrows():
        s=SeqRecord(Seq(r[seqkey]),id=r[idkey])
        seqs.append(s)
    return seqs

def dataframe_to_fasta(df, seqkey='translation', idkey='locus_tag',
                     descrkey='description',
                     outfile='out.faa'):
    """Genbank features to fasta file"""

    seqs=[]
    for i,row in df.iterrows():
        if descrkey in df.columns:
            d=row[descrkey]
        else:
            d=''
        rec = SeqRecord(Seq(row[seqkey]),id=row[idkey],
                            description=d)
        seqs.append(rec)
    SeqIO.write(seqs, outfile, "fasta")
    return outfile

def genbank_to_dataframe(infile, cds=False, quiet=True):
    """Get genome records from a genbank file into a dataframe
      returns a dataframe with a row for each cds/entry"""

    #genome = SeqIO.read(infile,'genbank')
    recs = list(SeqIO.parse(infile,'genbank'))
    genome = recs[0]
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
            if i in x:
                if type(x[i]) is list:
                    d[i] = x[i][0]
                else:
                    d[i] = x[i]
        allfeat.append(d)

    df = pd.DataFrame(allfeat,columns=featurekeys)
    df['length'] = df.translation.astype('str').str.len()
    #print (df)
    df = check_tags(df)
    if quiet == False:
        print('---- %s summary ----' %infile)
        s = genbank_summary(df)
        for i in s:
            print (i,':',s[i])
    if cds == True:
        df = get_cds(df)
        df['order'] = range(1,len(df)+1)
    #print (df)
    if len(df) == 0:
        print ('ERROR: genbank file return empty data, check that the file contains protein sequences '\
               'in the translation qualifier of each protein feature.' )
    return df

def check_tags(df):
    """Check genbank tags to make sure they are not empty"""

    def replace(x):
        if pd.isnull(x.locus_tag):
            return x.gene
        else:
            return x.locus_tag
    df['locus_tag'] = df.apply(replace,1)
    return df

def genbank_summary(df):
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

def find_keyword(f):
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

def index_genbank_features(gb_record, feature_type, qualifier):
    """Index features by qualifier value for easy access"""

    answer = dict()
    for (index, feature) in enumerate(gb_record.features):
        print (index, feature)
        if feature.type==feature_type:
            if qualifier in feature.qualifiers:
                values = feature.qualifiers[qualifier]
                if not type(values) is types.ListType:
                    values = [values]
                for value in values:
                    if value in answer:
                        print ("WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index))
                    else:
                        answer[value] = index
    return answer

def get_genes_by_location(genome, feature, within=20):
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

def get_translation(feature, genome, cds=True):
    """Check the translation of a cds feature"""

    trtable = "Bacterial"
    q = feature.qualifiers
    #trans = q1['translation'][0]
    seq = feature.extract(genome.seq)
    e=None
    try:
        protein = seq.translate(table=trtable,cds=cds,to_stop=True)
        #print ('protein seq:',protein)
    except Exception as e:
        protein = ''
    return protein, e

def pairwise_alignment(rec1,rec2):
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62
    gap_open = -10
    gap_extend = -0.5
    alns = pairwise2.align.localds(rec1, rec2, matrix, gap_open, gap_extend)
    return alns

def clustal_alignment(filename=None, seqs=None, command="clustalw"):
    """Align 2 sequences with clustal"""

    if filename == None:
        filename = 'temp.faa'
        SeqIO.write(seqs, filename, "fasta")
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import ClustalwCommandline
    cline = ClustalwCommandline(command, infile=filename)
    stdout, stderr = cline()
    align = AlignIO.read(name+'.aln', 'clustal')
    return align

def needle_alignment(seq1,seq2,outfile='needle.txt'):
    """Align 2 sequences with needle"""

    SeqIO.write(seq1, 'alpha.faa', "fasta")
    SeqIO.write(seq2, 'beta.faa', "fasta")
    from Bio.Emboss.Applications import NeedleCommandline
    cline = NeedleCommandline(asequence='alpha.faa', bsequence='beta.faa',
                              gapopen=30, gapextend=0.5, outfile=outfile)
    stdout, stderr = cline()
    align = AlignIO.read('needle.txt',"emboss")
    return align

def muscle_alignment(filename=None, seqs=None):
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

def show_alignment(aln, diff=False, offset=0):
    """
    Show a sequence alignment
        Args:
            aln: alignment
            diff: whether to show differences
    """

    ref = aln[0]
    l = len(aln[0])
    n=60
    chunks = [(i,i+n) for i in range(0, l, n)]
    for c in chunks:
        start,end = c
        lbls = np.arange(start,end,10)-offset
        print (('%-21s' %'name'),''.join([('%-10s' %i) for i in lbls]))
        print (('%21s' %ref.id[:20]), ref.seq[start:end])

        if diff == True:
            for a in aln[1:]:
                diff=''
                for i,j in zip(ref,a):
                    if i != j:
                        diff+=j
                    else:
                        diff+='-'
                name = a.id[:20]
                print (('%21s' %name), diff[start:end])
        else:
            for a in aln[1:]:
                name = a.id[:20]
                print (('%21s' %name), a.seq[start:end])
    return

def get_identity(aln):
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

def format_alignment(aln):
    t=''
    for i in range(0,len(aln[0]),80):
        for a in aln:
            t+=('%15s' %a.id) +' '+ a.seq.tostring()[i:i+80]+'\n'
        t+='\n'
    return t

def alignment_to_dataframe(aln):
    """Sequence alignment to dataframe"""

    alnrows = [[a.id,str(a.seq),a.description] for a in aln]
    df = pd.DataFrame(alnrows,columns=['name','seq','description'])
    return df

def get_feature_qualifier(f, qualifier):
    if f.qualifiers.has_key(qualifier):
        fq = f.qualifiers[qualifier][0].lower()
    else:
        fq = None
    return fq

def get_sequence(genome, name):
    """Get the sequence for a protein in a dataframe with
       genbank/sequence data"""
    return genome[genome.locus_tag==name].translation.iloc[0]

def fetch_protein_sequences(searchterm, filename='found.fa' ):
    """
    Fetch protein seqs using ncbi esearch and save results to a
    fasta file.
    Args:
        searchterm: entrez search term
        filename: fasta file name to save results
    Returns:
        sequence records as a dataframe
    """

    from Bio import Entrez
    from Bio import SeqIO
    Entrez.email = "A.N.Other@example.com"

    handle = Entrez.esearch(db="protein", term=searchterm, retmax=200)
    record = Entrez.read(handle)
    handle.close()
    #fetch the sequences
    handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=record['IdList'])
    seq_record = SeqIO.parse(handle, "fasta")
    recs = [r for r in seq_record]
    handle.close()
    outfile = open(filename, "w")
    SeqIO.write(recs, outfile, "fasta")

    df = fasta_to_dataframe(filename)
    #remove redundancy
    df = df.drop_duplicates('translation')
    df = df[-df.translation.str.contains('X')]
    print ('%s non-redundant sequences retrieved' %len(df))
    #save as fasta file
    dataframe_to_fasta(df, outfile=filename)
    return recs


def show_alignment_html(alnrows, seqs, width=80, fontsize=15, label='name'):
    """
    Get html display of sub-sequences on multiple protein alignment.
    Args:
        alnrows: a dataframe of aligned sequences
        seqs: sub-sequences/epitopes to draw if present
        label: key from dataframe to use as label for sequences
    Returns:
        html code
    """

    import matplotlib as mpl
    l=len(seqs[0])
    found = []
    for row in alnrows.seq:
        x = [row.find(s) for s in seqs]
        x = [i for i in x if i!=-1]
        #print x
        found.append(x)

    seqhtml=[]
    f=[]
    [f.extend(i) for i in found]
    f = sorted(list(set(f)))
    cmap = mpl.cm.get_cmap('Set3')
    c=1
    #unique color for each found sub-sequence
    colors={}
    for i in f:
        clr = cmap(float(c+0.1)/len(f))
        colors[i] = mpl.colors.rgb2hex(clr)
        c+=1

    seqhtml.append('<div style="font-family: monospace;letter-spacing: -.3em;font-size:%spx">' %fontsize)
    clr = ''
    chunks = []
    alnlength = len(alnrows.iloc[0].seq)
    l = 11
    for idx in range(0,alnlength,width):
        f=0
        seqhtml.append('<span style="letter-spacing:.2em;font-weight: bold">%s</span><br>' %idx)
        cidx=0
        for x,row in alnrows.iterrows():
            if len(found[f])==0:
                f+=1
                continue
            try:
                name = row[label]
            except:
                name = row.definition
            seq  = row.seq
            for i in range(idx,idx+width):
                if i>alnlength-1: continue
                if i in found[f]:
                    cidx = i
                    clr = colors[cidx]
                elif i-cidx >= l:
                    clr = ''
                seqhtml.append('<span style="background-color:%s">%s </span>' %(clr,seq[i]))

            clr = ''
            seqhtml.append('<span> &nbsp </span>')
            seqhtml.append('<span style="letter-spacing:.1em; font-weight: bold">%s </span>' %name)
            seqhtml.append('<br>')
            f+=1
    seqhtml = ''.join(seqhtml)
    return seqhtml