#!/usr/bin/env python

"""
    Utilities for docking/peptide stuff.
    Created March 2013
    Copyright (C) Damien Farrell
"""

import os, math, csv, string
import shutil
import Image, ImageFont, ImageDraw
import ConfigParser
import numpy as np
import pandas as pd
import matplotlib
import pylab as plt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import PDB

home = os.path.expanduser("~")

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

def addDicts(a, b):
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
        print src
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

def getSymmetricDataFrame(m):
    x = symmetrize(m)
    return pd.DataFrame(x, columns=m.columns,index=m.index)

'''def getcsvdata(f, delimiter=','):
    cr = csv.reader(open(f,'r'),delimiter=delimiter)
    a = [i for i in cr]
    return a'''

def parseConfig(conffile=None):
    """Parse the config file"""
    f = open(conffile,'r')
    cp = ConfigParser.ConfigParser()
    try:
        cp.read(conffile)
    except Exception,e:
        print 'failed to read config file! check format'
        print 'Error returned:', e
        return
    obj = setAttributesfromConfigParser(cp)
    return obj

def writeDefaultConfig(conffile='default.conf', defaults={}):
    """Write a default conf file"""
    if not os.path.exists(conffile):
        cp = createConfigParserfromDict(defaults, ['base'])
        cp.write(open(conffile,'w'))
    return conffile

def createConfigParserfromDict(data, sections, **kwargs):
    """Helper method to create a ConfigParser from a dict and/or keywords"""
    cp = ConfigParser.ConfigParser()
    for s in sections:
        cp.add_section(s)
        if not data.has_key(s):
            continue
        for i in data[s]:
            name,val = i
            cp.set(s, name, val)
    #use kwargs to create specific settings in the appropriate section
    for s in cp.sections():
        opts = cp.options(s)
        for k in kwargs:
            if k in opts:
                cp.set(s, k, kwargs[k])
    return cp

def setAttributesfromConfigParser(cp, obj=None):
    """A helper method that makes the options in a ConfigParser object
       attributes of an arbitrary object, obj """

    if obj == None:
        class Object(object):
            pass
        obj = Object()
    for s in cp.sections():
        obj.__dict__[s] = cp.items(s)
        for f in cp.items(s):
            try: val=int(f[1])
            except: val=f[1]
            obj.__dict__[f[0]] = val
    return obj

def getSequencefromFasta(filename):
    rec = SeqIO.parse(open(filename,'r'),'fasta').next()
    seq = str(rec.seq)
    return seq

def getListfromConfig(string, types='int'):
    """Extract a list from a comma separated config entry"""
    if types == 'int':
        vals = [int(i) for i in string.split(',')]
    return vals

def findFilefromString(files, string):
    for f in files:
        if string in os.path.splitext(f)[0]:
            return f
    return ''

def findFiles(path, ext='txt'):
    """List files in a dir of a specific type"""
    if not os.path.exists(path):
        print 'no such directory: %s' %path
        return []
    files=[]
    for dirname, dirnames, filenames in os.walk(path):
        for f in filenames:
            name = os.path.join(dirname, f)
            if f.endswith(ext):
                files.append(name)
    return files

def findFolders(path):
    if not os.path.exists(path):
        print 'no such directory: %s' %path
        return []
    dirs = []
    for dirname, dirnames, filenames in os.walk(path):
        dirs.append(dirname)
    return dirs

def reorderFilenames(files, order):
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

def combineImages(path=None, imgfiles=None, cols=3, size=300):
    """Combine several images in one on a grid"""

    font = ImageFont.truetype("Arial.ttf", 15)
    x=size
    w=20
    i=0; j=0
    if imgfiles == None:
        imgfiles = findFiles(path, 'png')
    width = cols*(x+w)
    height = int(math.ceil(float(len(imgfiles))/cols)*x)
    new_im = Image.new('RGBA', (width, height), 'white')
    for f in imgfiles:
        name = os.path.basename(f).split('.')[0]
        if not os.path.exists(f):
            continue
        im = Image.open(f)
        im.thumbnail((x,x))
        new_im.paste(im, (i*x+w,j*x+w))
        draw = ImageDraw.Draw(new_im)
        draw.text((i*x+w,j*x+w), name, (0,0,0), font=font)
        i+=1
        if i>=cols:
            i=0; j+=1
    #new_im.show()
    path = os.path.split(imgfiles[0])[0]
    new_im.save(os.path.join(path,"summary.png"))
    return

def removeNans(data):
    """Remove nans from lists"""
    for i in data[:]:
        ind = data.index(i)
        for j in i:
            if np.isnan(j):
                data.remove(i)
                break
    return data

def cluster(X=None, datalabels=None, nc=2):
    """Simple clustering for data of n-dimensions"""
    from sklearn.cluster import KMeans
    from sklearn.cluster import AffinityPropagation

    C = KMeans(n_clusters=nc,n_init=10,init='random')
    C.fit(X[:,:1])

    #C = AffinityPropagation(preference=-80,damping=0.5).fit(X)
    #cluster_centers_indices = C.cluster_centers_indices_

    clust = {}
    for (i, label) in enumerate(C.labels_):
        key = C.cluster_centers_[label][0]
        #print label,key, datalabels[i],X[i][1]
        if not clust.has_key(key):
            clust[key]=[]
        clust[key].append(datalabels[i])
    #print clust
    return C, clust

def plotCluster(X, est, ax=None):

    import pylab as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d import proj3d

    labels = est.labels_

    '''for x,y,n in zip(X[:, 0], X[:, 1], datalabels):
        x2, y2, _ = proj3d.proj_transform(x,y,1, ax.get_proj())
        if x>22:
            ax.annotate('{}'.format(n), xy=(x2,y2), xytext=(-5, 5), ha='left',
                    textcoords='offset points', fontsize=8)'''
    if ax == None:
        fig = plt.figure()
        ax=fig.add_subplot(111,projection='3d')
        ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=labels.astype(np.float), alpha=.6)
    else:
        ax.scatter(X[:, 0], X[:, 1], c=labels.astype(np.float), s=30,alpha=.6)
        #ax.set_xlabel('mean rmsd')
        #ax.set_ylabel('binding energy(kcal/mol)')
    #fig.savefig('cluster.png')
    return ax

def plotSurface(X):
    """Plot surface for 3 sets of arrays"""
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d import proj3d
    f=plt.figure()
    ax=f.add_subplot(111,projection='3d')
    xi=np.arange(10,14,0.05)
    yi=np.arange(12,16,0.05)
    z = matplotlib.mlab.griddata(X[:,0], X[:,1], X[:,2], xi, yi, interp='nn')
    x, y = np.meshgrid(xi, yi)
    ax.plot_surface(x, y, z)
    return f

def drawBoxes():
    from matplotlib.patches import Rectangle
    f=plt.figure()
    ax=f.add_subplot(111)
    epos = []
    seq=range(len(seqstr))
    y=np.random.randint(1,10,len(seqstr))
    ax.plot(seq,y)
    for e in exp:
        seq=e[0]
        i = seqstr.find(seq)
        print e,i
        #epos.append(i)
        if e[1] > 0.5:
            polygon = Rectangle((i,2), 1,0.5, color='r')
            #patches.append(polygon)
            ax.add_patch(polygon)
    f.savefig('testboxes.png')

def readIEDB(filename, key='Epitope ID'):
    """Load iedb peptidic csv file and return dataframe"""
    #cr = csv.reader(open(filename,'r'))
    cr = csv.DictReader(open(filename,'r'),quotechar='"')
    cr.fieldnames = [field.strip() for field in cr.fieldnames]
    D={}
    for r in cr:
        k = r[key]
        D[k] = r
    return D

def getCSVData(filename, key):
    """Get data from cvs file into dict, also returns
       original order of key names in file as a list"""
    cr = csv.DictReader(open(filename,'r'))
    data = {}
    order = []
    for r in cr:
        k = r[key]
        data[k] = r
        order.append(k)
    fields = cr.fieldnames
    return data, order

def fetchPDB(name, path):
    """Fetch a pdb and save to path"""
    from Bio.PDB import PDBList
    pdbname = os.path.join(path,name+'.pdb')
    pdbl = PDBList()
    filename = pdbl.retrieve_pdb_file(name,pdir=path)
    os.rename(filename, pdbname)
    return

def getSequencefromPDB(pdbfile, chain='C', index=0):
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

def renumberChain(chain):
    """Renumber pdb chain from 1, if it starts from another number"""
    first = chain.child_list[0].id[1]
    #print chain, first
    if first == 0:
        offset = -1
    elif first > 1:
        offset = first-1
    else:
        offset = 0
    for residue in chain.child_list:
        i=residue.id[1]
        residue.id = (' ', i-offset, ' ')
    return

def resetChainNumbering(chain):
    """Renumber from 1, regardless of previous numbering"""
    i=1
    for residue in chain.child_list[:]:
        residue.id = (' ', i, ' ')
        i+=1
    return

def saveStructure(structure, filename):
    w = PDB.PDBIO()
    w.set_structure(structure)
    w.save(filename, write_end=1)
    return filename

def removeChains(pdbfile, start=3):
    print pdbfile
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdbfile, pdbfile)
    model = structure[0]
    remove = string.uppercase[3:]
    for chain in model.child_list[:]:
        if chain.id in remove:
            model.detach_child(chain.id)
        else:
            for residue in chain.child_list[:]:
                if residue.id[0] != " ":
                    chain.detach_child(residue.id)
    renumberChain(model['C'])
    return structure

def preparePDB(pdbfile, start=3):
    """Prepare an MHC pdb for use in docking by removing unused atoms,
       combining A+B chains and renumbering chain C"""

    print pdbfile
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdbfile, pdbfile)
    model = structure[0]
    remove = string.uppercase[3:]
    for chain in model.child_list[:]:
        if chain.id in remove:
            model.detach_child(chain.id)
        else:
            for residue in chain.child_list[:]:
                #remove all non-standard residues
                if residue.id[0] != " ":
                    chain.detach_child(residue.id)

    removeresidues = [('A', range(86,200)), ('B', range(96,200))]
    for c in removeresidues:
        chain = model[c[0]]
        remove = c[1]
        for residue in chain.child_list[:]:
            id = residue.id
            if id[1] in remove:
                chain.detach_child(id)
    #renumber chain A
    chaina = model['A']
    renumberChain(chaina)
    #renumber chain B
    chainb = model['B']
    i = chaina.child_list[-1].id[1] + 1
    for residue in chainb:
        residue.id = (' ', i, ' ')
        i+=1
    chainb.id = 'A'
    #model.detach_child('B')
    #print '%s chains, length chain a=%s' %(len(model.child_list), len(chaina.child_list))
    #renumber chain c
    renumberChain(model['C'])
    w = PDB.PDBIO()
    w.set_structure(structure)
    name = os.path.splitext(pdbfile)[0]
    filename = name+'_alt.pdb'
    w.save(filename, write_end=1)
    #print 'saved file %s' %filename

    #save unbound
    model.detach_child('C')
    filename2 = name+'_unbound.pdb'
    w.save(filename2, write_end=1)
    return filename, filename2

def prepareNmer(pdbfile, length=9, start=3, chain='C'):
    """Cut and renumber the n-mer peptide from the native complex"""

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdbfile, pdbfile)
    model = structure[0]
    ch = model[chain]
    end = start+length
    for residue in ch.child_list[:]:
        pos = residue.id[1]
        if pos < start or pos >= end:
            ch.detach_child(residue.id)
    renumberChain(ch)
    w = PDB.PDBIO()
    w.set_structure(structure)
    name = os.path.splitext(pdbfile)[0]
    filename = name+'%smer.pdb' %length
    w.save(filename)
    print 'saved n-mer file %s' %filename
    return filename

def miscFixes(pdbfile):
    """Misc fixes to pdb"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdbfile, pdbfile)
    model = structure[0]
    chc = model['C']
    resetChainNumbering(chc)
    w = PDB.PDBIO()
    w.set_structure(structure)
    w.save(pdbfile)
    return

def cleanHaddockPDB(pdbfile, new=None):
    """Removes trailing chain id from haddock pdbs"""
    f = open(pdbfile, 'r')
    if new == None:
        new = pdbfile
    lines = list(f)
    f.close()
    n = open(new,'w')
    for line in lines:
        if line.startswith('ATOM'):
            n.write(line[:-7]+'\n')
    return

def filterIEDBFile(filename, field, search):
    """Return filtered iedb data"""
    X = pd.read_csv(filename)
    cols = ['PubMed ID','Author','Journal','Year','T Cell ID','MHC Allele Name',
                'Epitope Linear Sequence','Epitope Source Organism Name']
    y = X[X[field].str.contains(search)]
    print y[cols]
    y.to_csv('filtered.csv',cols=cols)
    return y

def test():
    sourcefasta = os.path.join(home,'dockingdata/fastafiles/1KLU.fasta')
    findClosestStructures(sourcefasta)
    #fetchPDBList('MHCII_homologs.csv')

if __name__ == '__main__':
    test()
