#!/usr/bin/env python

"""
    Utilities for epitopepredict
    Created March 2013
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import os, math, csv, string
import shutil
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import PDB

home = os.path.expanduser("~")

def venndiagram(names, labels, ax=None, colors=('r','g','b')):
    """Plot a venn diagram"""

    from matplotlib_venn import venn2,venn3
    import pylab as plt
    f=None
    if ax==None:
        f=plt.figure(figsize=(4,4))
        ax=f.add_subplot(111)
    if len(names)==2:
        n1,n2=names
        v = venn2([set(n1), set(n2)], set_labels=labels, set_colors=colors)
    elif len(names)==3:
        n1,n2,n3=names
        v = venn3([set(n1), set(n2), set(n3)], set_labels=labels, set_colors=colors)
    ax.axis('off')
    #f.patch.set_visible(False)
    ax.set_axis_off()
    return f

def compress(filename, remove=False):
    """Compress a file with gzip"""
    import gzip
    fin = open(filename, 'rb')
    fout = gzip.open(filename+'.gz', 'wb')
    fout.writelines(fin)
    fout.close()
    fin.close()
    if remove == True:
        os.remove(filename)
    return

def rmse(ar1, ar2):
    """Mean squared error"""
    ar1 = np.asarray(ar1)
    ar2 = np.asarray(ar2)
    dif = ar1 - ar2
    dif *= dif
    return np.sqrt(dif.sum()/len(ar1))

def add_dicts(a, b):
    return dict((n, a.get(n, 0)+b.get(n, 0)) for n in set(a)|set(b))

def copyfile(source, dest, newname=None):
    """Helper method to copy files"""

    if not os.path.exists(source):
        #print 'no such file %s' %source
        return False
    shutil.copy(source, newname)
    dest = os.path.join(dest, newname)
    if os.path.exists(dest):
        os.remove(dest)
    shutil.move(newname, dest)
    return True

def copyfiles(path, files):
    for f in files:
        src = os.path.join(path, f)
        print (src)
        if not os.path.exists(src):
            return False
        shutil.copy(src, f)
    return True

def symmetrize(m, lower=True):
    """Return symmetric array"""
    m=m.fillna(0)
    if lower==True:
        return np.tril(m) + np.triu(m.T) - np.diag(np.diag(m))
    else:
        return np.triu(m) + np.tril(m.T) - np.diag(np.diag(m))

def get_symmetric_data_frame(m):
    x = symmetrize(m)
    return pd.DataFrame(x, columns=m.columns,index=m.index)

def find_filefrom_string(files, string):
    for f in files:
        if string in os.path.splitext(f)[0]:
            return f
    return ''

def find_files(path, ext='txt'):
    """List files in a dir of a specific type"""
    if not os.path.exists(path):
        print ('no such directory: %s' %path)
        return []
    files=[]
    for dirname, dirnames, filenames in os.walk(path):
        for f in filenames:
            name = os.path.join(dirname, f)
            if f.endswith(ext):
                files.append(name)
    return files

def find_folders(path):
    if not os.path.exists(path):
        print ('no such directory: %s' %path)
        return []
    dirs = []
    for dirname, dirnames, filenames in os.walk(path):
        dirs.append(dirname)
    return dirs

def reorder_filenames(files, order):
    """reorder filenames by another list order(seqs)"""
    new = []
    for i in order:
        found=False
        for f in files:
            if i in f:
                new.append(f)
                found=True
        if found==False:
            new.append('')
    return new

def read_iedb(filename, key='Epitope ID'):
    """Load iedb peptidic csv file and return dataframe"""
    #cr = csv.reader(open(filename,'r'))
    cr = csv.DictReader(open(filename,'r'),quotechar='"')
    cr.fieldnames = [field.strip() for field in cr.fieldnames]
    D={}
    for r in cr:
        k = r[key]
        D[k] = r
    return D

def get_sequencefrom_pdb(pdbfile, chain='C', index=0):
    """Get AA sequence from PDB"""
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure(pdbfile,pdbfile)
    ppb = PDB.PPBuilder()
    model = struct[0]
    peptides = ppb.build_peptides(model[chain])
    seq=''
    for i,pep in enumerate(peptides):
        seq+=str(pep.get_sequence())
    return seq

def filter_iedb_file(filename, field, search):
    """Return filtered iedb data"""
    X = pd.read_csv(filename)
    cols = ['PubMed ID','Author','Journal','Year','T Cell ID','MHC Allele Name',
                'Epitope Linear Sequence','Epitope Source Organism Name']
    y = X[X[field].str.contains(search)]
    print (y[cols])
    y.to_csv('filtered.csv',cols=cols)
    return y

def search_pubmed(term, max_count=100):

    from Bio import Entrez
    from Bio import Medline

    def fetch_details(id_list):
        ids = ','.join(id_list)
        Entrez.email = 'your.email@example.com'
        handle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                               id=ids)
        results = Entrez.read(handle)
        return results

    def search(query):
        Entrez.email = 'your.email@example.com'
        handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax=max_count,
                                retmode='xml',
                                term=query)
        results = Entrez.read(handle)
        return results

    results = search(term)
    id_list = results['IdList']
    papers = fetch_details(id_list)
    for i, paper in enumerate(papers):
        print("%d) %s" % (i+1, paper['MedlineCitation']['Article']['ArticleTitle']))
        # Pretty print the first paper in full to observe its structure
        #import json
        #print(json.dumps(papers[0], indent=2, separators=(',', ':')))

def test():
    sourcefasta = os.path.join(home,'dockingdata/fastafiles/1KLU.fasta')
    findClosestStructures(sourcefasta)
    #fetchPDBList('MHCII_homologs.csv')

if __name__ == '__main__':
    test()
