#!/usr/bin/env python

"""
    MHCpredict plotting
    Created February 2016
    Copyright (C) Damien Farrell
"""

import sys, os
from collections import OrderedDict
import pylab as plt
import numpy as np
import pandas as pd
from . import base

colormaps={'tepitope':'Greens','netmhciipan':'Oranges','iedbmhc2':'Pinks',
               'threading':'Purples','iedbmhc1':'Blues'}
colors = {'tepitope':'green','netmhciipan':'orange',
           'iedbmhc1':'blue','iedbmhc2':'pink','threading':'purple'}

def plotTracks(preds, title='', n=2, width=820, height=None,
                seqdepot=None, bcell=None, exp=None, tools=True):
    """Plot binding predictions as parallel tracks of blocks for each allele.
       This uses Bokeh.
       returns: a bokeh figure for embedding or displaying in a notebook"""

    from collections import OrderedDict
    from bokeh.models import Range1d,HoverTool,FactorRange,Grid,GridPlot,ColumnDataSource
    from bokeh.plotting import Figure
    import matplotlib as mpl

    if tools == True:
        tools="xpan, xwheel_zoom, resize, hover, reset, save"
    else:
        tools=''

    alls=1
    for m in preds:
        alls += len(preds[m].data.groupby('allele'))
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

    #lists for hover data
    #we plot all rects at once
    x=[];y=[];allele=[];widths=[];clrs=[];peptide=[]
    predictor=[];position=[];score=[];leg=[]
    l=80
    for m in preds:
        pred = preds[m]
        cmap = mpl.cm.get_cmap(colormaps[m])
        df = pred.data
        sckey = pred.scorekey
        pb = pred.getPromiscuousBinders(data=df,n=n)
        if len(pb) == 0:
            continue
        l = pred.getLength()
        grps = df.groupby('allele')
        alleles = grps.groups.keys()
        if len(pb)==0:
            continue
        c=colors[m]
        leg.append(m)

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

    seqlen = pred.data.pos.max()+l
    plot.set(x_range=Range1d(start=0, end=seqlen+1))
    plot.xaxis.major_label_text_font_size = "8pt"
    plot.xaxis.major_label_text_font_style = "bold"
    plot.ygrid.grid_line_color = None
    plot.yaxis.major_label_text_font_size = '0pt'
    plot.xaxis.major_label_orientation = np.pi/4
    return plot
