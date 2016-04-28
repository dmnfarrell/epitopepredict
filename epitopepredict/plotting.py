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
colors = {'tepitope':'green','netmhciipan':'orange',
           'iedbmhc1':'blue','iedbmhc2':'pink','threading':'purple'}

def plot_tracks(preds, title='', n=2, cutoff_method='default', name=None,
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

def mpl_plot_tracks(preds, name, cldist=7, n=2, cutoff_method='default',
                legend=False, figsize=(13,4), ax=None):
    """Plot binders as bars per allele - defunct"""

    from matplotlib.patches import Rectangle
    if ax==None:
        fig=plt.figure(figsize=figsize)
        ax=fig.add_subplot(111)

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
        pb = pred.getPromiscuousBinders(data=df, n=n,
                                        cutoff_method=cutoff_method)
        if len(pb) == 0:
            continue
        l = base.getLength(pb)
        seqlen = df.pos.max()+l
        grps = df.groupby('allele')
        alleles.extend(grps.groups.keys())
        if len(pb)==0:
            continue
        c=colors[m]
        leg.append(m)

        for a,g in grps:
            b = pred.getBinders(data=g)
            b = b[b.pos.isin(pb.pos)] #only promiscuous
            b.sort_values('pos',inplace=True)
            #scores = b[sckey].values
            pos = b['pos'].values
            for x in pos:
                rect = ax.add_patch(Rectangle((x,y), l, 1, facecolor=c, lw=1.5, alpha=0.7))
            y+=1
        handles.append(rect)

    ax.set_xlim(0, seqlen)
    ax.set_ylim(0, len(alleles))
    ax.set_ylabel('allele')
    ax.set_yticks(np.arange(.5,len(alleles)+.5))
    ax.set_yticklabels(alleles)
    ax.grid(b=True, which='major', alpha=0.5)
    ax.set_title(name, fontsize=15)
    if legend == True:
        ax.legend(handles, leg, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=3)
    #plt.tight_layout()
    return ax

def mpl_plot_regions(coords, ax, color='red', label=''):
    """Highlight regions in a prot binder plot"""

    from matplotlib.patches import Rectangle
    #l = len(seqs.head(1)['key'].max())
    h = ax.get_ylim()[1]
    for c in coords:
        x,l = c
        ax.add_patch(Rectangle((x,0), l, h,
            facecolor=color, lw=.5, alpha=0.5))
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
