#!/usr/bin/env python

"""
    MHCpredict plotting
    Created February 2016
    Copyright (C) Damien Farrell
"""

import sys, os
from collections import OrderedDict
import pylab as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from . import base

colormaps={'tepitope':'Greens','netmhciipan':'Oranges','iedbmhc2':'Pinks',
               'threading':'Purples','iedbmhc1':'Blues'}
defaultcolors = {'tepitope':'green','netmhciipan':'orange',
           'iedbmhc1':'blue','iedbmhc2':'pink','threading':'purple'}

def plot_heatmap(df, ax=None, figsize=(6,6), **kwargs):
    """Plot a heatmap - move to plotting?"""

    if ax==None:
        fig=plt.figure(figsize=figsize)
        ax=fig.add_subplot(111)
    else:
        fig = ax.get_figure()
    df = df._get_numeric_data()
    hm = ax.pcolor(df, **kwargs)
    #fig.colorbar(hm, ax=ax)
    ax.set_xticks(np.arange(0.5, len(df.columns)))
    ax.set_yticks(np.arange(0.5, len(df.index)))
    ax.set_xticklabels(df.columns, minor=False, fontsize=14,rotation=90)
    ax.set_yticklabels(df.index, minor=False, fontsize=14)
    ax.set_ylim(0, len(df.index))
    hm.set_clim(0,1)
    #ax.grid(True)
    plt.tight_layout()
    return ax

def bokeh_plot_tracks(preds, title='', n=2, cutoff_method='default', name=None,
                width=820, height=None, tools=True,
                seqdepot=None, bcell=None, exp=None):
    """
    Plot binding predictions as parallel tracks of blocks for each allele.
    This uses Bokeh.
    Args:
        title: plot title
        n: min alleles to display
        name: name of protein to show if more than one in data

    Returns: a bokeh figure for embedding or displaying in a notebook
    """

    from collections import OrderedDict
    from bokeh.models import Range1d,HoverTool,FactorRange,Grid,GridPlot,ColumnDataSource
    from bokeh.plotting import Figure
    import matplotlib as mpl

    if tools == True:
        tools="xpan, xwheel_zoom, resize, hover, reset, save"
    else:
        tools=''

    alls=1
    for p in preds:
        if p.data is None:
            continue
        alls += len(p.data.groupby('allele'))
    if height==None:
        height = 130+10*alls
    yrange = Range1d(start=0, end=alls+3)
    plot = Figure(title=title,title_text_font_size="11pt",plot_width=width,
                  plot_height=height, y_range=yrange,
                y_axis_label='allele',
                tools=tools,
                background_fill="#FAFAFA",
                toolbar_location="below")
    h=3
    if bcell != None:
        plotBCell(plot, bcell, alls)
    if seqdepot != None:
        plotAnnotations(plot,seqdepot)
    if exp is not None:
        plotExp(plot, exp)

    #plotRegions(plot)

    #we plot all rects at once
    x=[];y=[];allele=[];widths=[];clrs=[];peptide=[]
    predictor=[];position=[];score=[];leg=[]
    l=80
    for pred in preds:
        #pred = preds[m]
        m = pred.name
        cmap = mpl.cm.get_cmap(colormaps[m])
        df = pred.data
        if df is None or len(df) == 0:
            print('no data to plot for %s' %m)
            continue
        if name != None:
            df = df[df.name==name]
        sckey = pred.scorekey
        pb = pred.getPromiscuousBinders(data=df,n=n, cutoff_method=cutoff_method)
        if len(pb) == 0:
            continue
        l = base.getLength(pb)
        grps = df.groupby('allele')
        alleles = grps.groups.keys()
        if len(pb)==0:
            continue
        c=colors[m]
        leg.append(m)
        seqlen = df.pos.max()+l

        for a,g in grps:
            b = pred.getBinders(data=g)
            b = b[b.pos.isin(pb.pos)] #only promiscuous
            b.sort_values('pos',inplace=True)
            scores = b[sckey].values
            score.extend(scores)
            pos = b['pos'].values
            position.extend(pos)
            x.extend(pos+(l/2.0)) #offset as coords are rect centers
            widths.extend([l for i in scores])
            clrs.extend([c for i in scores])
            y.extend([h+0.5 for i in scores])
            alls = [a for i in scores]
            allele.extend(alls)
            peptide.extend(list(b.peptide.values))
            predictor.extend([m for i in scores])
            h+=1

    source = ColumnDataSource(data=dict(x=x,y=y,allele=allele,peptide=peptide,
                                    predictor=predictor,position=position,score=score))
    plot.rect(x,y, width=widths, height=0.8,
         #x_range=Range1d(start=1, end=seqlen+l),
         color=clrs,line_color='gray',alpha=0.7,source=source)

    hover = plot.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("allele", "@allele"),
        ("position", "@position"),
        ("peptide", "@peptide"),
        ("score", "@score"),
        ("predictor", "@predictor"),
    ])

    plot.set(x_range=Range1d(start=0, end=seqlen+1))
    plot.xaxis.major_label_text_font_size = "8pt"
    plot.xaxis.major_label_text_font_style = "bold"
    plot.ygrid.grid_line_color = None
    plot.yaxis.major_label_text_font_size = '0pt'
    plot.xaxis.major_label_orientation = np.pi/4
    return plot

def plot_tracks(preds, name, n=2, cutoff=5, value='score',
                legend=False, colormap='Paired', figsize=None, ax=None):
    """
    Plot binders as bars per allele using matplotlib.
    Args:
        preds: list of one or more predictors
        name: name of protein to plot
        n: number of alleles binder should be found in to be displayed
        cutoff: percentile cutoff to determine binders to show

    """

    import matplotlib as mpl
    from matplotlib.patches import Rectangle
    if ax==None:
        if figsize==None:
            h = sum([len(p.data.groupby('allele')) for p in preds])
            w=10
            h = round(h*.1+2)
            figsize = (w,h)
        fig=plt.figure(figsize=figsize)
        ax=fig.add_subplot(111)

    p = len(preds)
    cmap = mpl.cm.get_cmap(colormap)
    colors = { preds[i].name : cmap(float(i)/p) for i in range(p) }

    alleles = []
    leg = []
    y=0
    handles = []
    for pred in preds:
        m = pred.name
        df = pred.data
        if df is None or len(df) == 0:
            print('no data to plot for %s' %m)
            continue
        if name != None:
            df = df[df.name==name]
        sckey = pred.scorekey

        binders = pred.getBinders(name, cutoff=cutoff, value=value)
        #pass binders so it's not recalculated
        pb = pred.promiscuousBinders(binders=binders, n=n, value=value)

        if len(pb) == 0:
            continue
        l = base.getLength(pb)
        seqlen = df.pos.max()+l
        #print (m,df.pos.max())
        grps = df.groupby('allele')
        if m in colors:
            c=colors[m]
        else:
            c='blue'
        leg.append(m)
        order = sorted(grps.groups)
        alleles.extend(order)
        #for a,g in grps:
        for a in order:
            g = grps.groups[a]
            b = binders[binders.allele==a]
            b = b[b.pos.isin(pb.pos)] #only promiscuous
            b.sort_values('pos',inplace=True)
            pos = b['pos'].values+1 #assumes pos is zero indexed
            #clrs = [scmap.to_rgba(i) for i in b[sckey]]
            #for x,c in zip(pos,clrs):
            for x in pos:
                rect = ax.add_patch(Rectangle((x,y), l, 1, facecolor=c,
                                    lw=1.5, alpha=0.6))
            y+=1
        handles.append(rect)
    if len(leg) == 0:
        return
    ax.set_xlim(0, seqlen)
    ax.set_ylim(0, len(alleles))
    w=20
    if seqlen>500: w=100
    ax.set_xticks(np.arange(0, seqlen, w))
    ax.set_ylabel('allele')
    ax.set_yticks(np.arange(.5,len(alleles)+.5))
    fsize = 14-1*len(alleles)/40.
    ax.set_yticklabels(alleles, fontsize=fsize )
    ax.grid(b=True, which='major', alpha=0.5)
    ax.set_title(name, fontsize=16, loc='right')
    if legend == True:
        ax.legend(handles, leg, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=3)
    #plt.tight_layout()
    return ax

def mpl_plot_regions(coords, ax, color='red', label='', alpha=0.6):
    """Highlight regions in a prot binder plot"""

    from matplotlib.patches import Rectangle
    #l = len(seqs.head(1)['key'].max())
    h = ax.get_ylim()[1]
    for c in coords:
        x,l = c
        ax.add_patch(Rectangle((x,0), l, h,
            facecolor=color, lw=.8, alpha=alpha, zorder=0))
    return

def mpl_draw_labels(labels, coords, ax):
    """Add labels on axis"""

    bbox_args = dict(boxstyle='square',fc='whitesmoke')
    from matplotlib.transforms import blended_transform_factory
    tform = blended_transform_factory(ax.transData, ax.transAxes)
    for text, x in zip(labels,coords):
        xy = (x,-.05)
        an = ax.annotate(text, xy=xy, xycoords=tform,
                   ha='center', va="center",
                   size=14,
                   zorder=10, textcoords='offset points',
                   bbox=bbox_args)
    plt.subplots_adjust(bottom=0.1)
    return

def plot_bars(P, name, chunks=1, how='median', perc=0.8, color='black'):
    """
    Bar plots for sequence using median/mean/total scores.
    Args:
        P: predictor with data
        name: name of protein sequence
        chunks: break sequence up into 1 or more chunks
        how: method to calculate score bar value
        perc: percentile cutoff to show peptide
    """

    import seaborn as sns
    df = P.data[P.data.name==name].sort_values('pos')
    w = 10
    l= base.getLength(df)
    seqlen = df.pos.max()+l
    f,axs = plt.subplots(chunks,1,figsize=(15,2+2.5*chunks))
    if chunks == 1:
        axs = [axs]
    else:
        axs = list(axs.flat)

    funcs = {'median': np.median, 'mean': np.mean, 'sum': np.sum}
    #cb = P.consensusRankedBinders()
    #cb = cb[cb['rank']<20]

    grps = df.groupby('pos')
    key = P.scorekey
    X = grps.agg({key: np.median, 'peptide': base.first})
    cutoff = X[key].quantile(perc)
    X[key][X[key]<cutoff] = np.nan
    seqlist = X.peptide.apply( lambda x : x[0])
    seqchunks = np.array_split(X.index, chunks)

    for c in range(chunks):
        ax = axs[c]
        st = seqchunks[c][0]
        end = seqchunks[c][-1]
        p = X[st:end]
        #p = p[p.peptide.isin(cb.peptide)]
        ax.bar(p.index, p[key], width=1, color=color)
        ax.set_title(str(st)+'-'+str(end), loc='right')
        xseq = seqlist[st:end]
        if len(xseq)<150:
            fsize = 16-1*len(xseq)/20.
            ax.set_xlim(st,end)
            ax.set_xticks(p.index+0.5)
            ax.set_xticklabels(xseq, rotation=0, fontsize=fsize)
        ax.set_ylim(X[key].min(), X[key].max())
    f.suptitle(name+' - '+P.name)
    plt.tight_layout()
    return axs

def mpl_plot_seqdepot(annotation, ax):
    """Plot sedepot annotations"""

    from matplotlib.patches import Rectangle
    y=-1.5
    fontsize=12

    if 'signalp' in annotation:
        bbox_args = dict(boxstyle='rarrow', fc='white', lw=1, alpha=0.8)
        pos = annotation['signalp'].values()
        print (pos)
        for x in pos:
            an = ax.annotate('SP', xy=(x,y), xycoords='data',
                        ha='left', va="center", bbox=bbox_args,
                        size=fontsize)
    if 'tmhmm' in annotation:
        vals = annotation['tmhmm']
        pos = [i[0]+(i[1]-i[0])/2.0 for i in vals]
        widths = [i[1]-i[0] for i in vals]
        bbox_args = dict(boxstyle='round', fc='deepskyblue', lw=1, alpha=0.8)
        for x,w in zip(pos,widths):
            an = ax.annotate('TMHMM', xy=(x,y), xycoords='data',
                       ha='left', va="center", bbox=bbox_args,
                       size=fontsize)
    if 'pfam27' in annotation:
        vals = annotation['pfam27']
        text = [i[0] for i in vals]
        pos = [i[1]+(i[2]-i[1])/2.0 for i in vals]
        widths = [i[2]-i[1] for i in vals]
        print (pos,widths,text)
        bbox_args = dict(boxstyle='round', fc='coral', lw=1, alpha=0.8)
        for x,w,t in zip(pos,widths,text):
            an = ax.annotate(t, xy=(x,y), xycoords='data',
                       ha='left', va="center", bbox=bbox_args,
                       size=fontsize)

    ax.set_ylim(y-1, ax.get_ylim()[1])
    return

def plot_multiple(preds, names, kind='tracks', regions=None, genome=None, **kwargs):
    """Plot results for multiple proteins"""

    for prot in names:
        if kind == 'tracks':
            ax = plot_tracks(preds,name=prot,**kwargs)
        elif kind == 'bar':
            axs = plot_bars(preds[0],name=prot)
            ax = axs[0]
        if regions is not None:
            r = regions[regions.name==prot]
            print (r)
            #print genome[genome.locus_tag==prot]
            coords = (list(r.start),list(r.end-r.start))
            coords = zip(*coords)
            mpl_plot_regions(coords, ax, color='gray')
        #labels = list(r.peptide)
        #plotting.mpl_draw_labels(labels, coords, ax)
        if genome is not None:
            p = genome[genome['locus_tag']==prot]
            seq = p.translation.iloc[0]
            sd = analysis.getSeqDepot(seq)['t']
            mpl_plot_seqdepot(sd, ax)
        plt.tight_layout()
        plt.show()
    return

def plot_binder_map(P, name, values='rank', cutoff=20, chunks=1, cmap=None):
    """
    Plot heatmap of binders above a cutoff by rank or score.
    Args:
        P: predictor object with data
        name: name of protein to plot
        values: data column to use for plot data, 'score' or 'rank'
        cutoff: cutoff if using rank as values
        chunks: number of plots to split the sequence into
    """

    import seaborn as sns
    df = P.data[P.data.name==name].sort('pos')
    w = 10
    l= base.getLength(df)
    seqlen = df.pos.max()+l
    f,axs = plt.subplots(chunks,1,figsize=(15,3+2.5*chunks))
    if chunks == 1:
        axs = [axs]
    else:
        axs = list(axs.flat)

    if values == 'score':
        values = P.scorekey
        if cmap == None: cmap='RdBu_r'
    X = df.pivot_table(index='allele', columns='pos', values=values)
    if values == P.scorekey:
        #normalise score across alleles for clarity
        zscore = lambda x: (x - x.mean()) / x.std()
        X = X.apply(zscore, 1)
    if values == 'rank':
        X[X > cutoff] = 0
        if cmap == None: cmap='Blues'
    s = df.drop_duplicates(['peptide','pos'])
    seqlist = s.set_index('pos').peptide.apply( lambda x : x[0])
    #print seqlist
    seqchunks = np.array_split(X.columns, chunks)
    for c in range(chunks):
        ax = axs[c]
        p = X[seqchunks[c]]
        #plot heatmap
        sns.heatmap(p, cmap=cmap, cbar_kws={"shrink": .5},
                    vmin=min(X.min()), vmax=max(X.max()),
                    ax=ax, xticklabels=20)
        #show sequence on x-axis
        st = seqchunks[c][0]
        end = seqchunks[c][-1]
        xseq = seqlist[st:end]
        ax.set_title(str(st)+'-'+str(end), loc='right')
        ax.spines['bottom'].set_visible(True)
        if len(seqchunks[c])<150:
            fsize = 16-1*len(seqchunks[c])/20.
            ax.set_xticks(np.arange(0,len(xseq))+0.5)
            ax.set_xticklabels(xseq, rotation=0, fontsize=fsize)

    f.suptitle(name+' - '+P.name)
    plt.tight_layout()
    return ax
