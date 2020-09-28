#!/usr/bin/env python

"""
    epitopepredict plotting
    Created February 2016
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
import sys, os, math
from collections import OrderedDict
try:
    import matplotlib
    matplotlib.use('agg', warn=False)
    import pylab as plt
except:
    pass
import numpy as np
import pandas as pd
from . import base

colormaps={'tepitope':'Greens','netmhciipan':'Oranges','iedbmhc2':'Pinks',
            'iedbmhc1':'Blues'}
defaultcolors = {'tepitope':'green','netmhciipan':'orange','basicmhc1':'yellow',
                 'iedbmhc1':'blue','iedbmhc2':'pink'}


def plot_heatmap(df, ax=None, figsize=(6,6), **kwargs):
    """Plot a generic heatmap """

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

def get_seq_from_binders(P, name=None):
    """Get sequence from binder data. Probably better to store the sequences in the object?"""

    if P.data is None or len(P.data)==0:
        return
    if name is not None:
        data=P.data[P.data.name==name]
    else:
        data=P.data
    l = len(data.iloc[0].peptide)
    seqlen = data.pos.max()+l
    return seqlen

def get_bokeh_colors(palette='Set1'):

    from bokeh.palettes import brewer
    n = len(base.predictors)
    pal = brewer[palette][n]
    i=0
    clrs={}
    for m in base.predictors:
        clrs[m] = pal[i]
        i+=1
    return clrs

def get_sequence_colors(seq):
    """Get colors for a sequence"""

    from bokeh.palettes import brewer, viridis, plasma
    from Bio.PDB.Polypeptide import aa1
    pal = plasma(20)
    pal.append('white')
    aa1 = list(aa1)
    aa1.append('-')
    pcolors = {i:j for i,j in zip(aa1,pal)}
    text = list(seq)
    clrs =  {'A':'red','T':'green','G':'orange','C':'blue','-':'white'}
    try:
        colors = [clrs[i] for i in text]
    except:
        colors = [pcolors[i] for i in text]
    return colors

def bokeh_test(n=20,height=400):

    from bokeh.models import ColumnDataSource
    from bokeh.plotting import figure
    from bokeh.models.glyphs import Text, Rect, Circle
    data = {'x_values': np.random.random(n),
            'y_values': np.random.random(n)}
    source = ColumnDataSource(data=data)
    tools = "pan,wheel_zoom,hover,tap,reset,save"
    p = figure(plot_height=height,tools=tools)
    c = Circle(x='x_values', y='y_values', radius=.02, line_color='black', fill_color='blue', fill_alpha=.6)
    #p.circle(x='x_values', y='y_values', radius=.02, line_color='black', fill_color='blue', fill_alpha=.6, source=source)
    p.add_glyph(source, c)
    return p

def bokeh_summary_plot(df, savepath=None):
    """Summary plot"""

    from bokeh.plotting import figure
    from bokeh.layouts import column
    from bokeh.models import ColumnDataSource,Range1d,HoverTool,TapTool,CustomJS,OpenURL

    TOOLS = "pan,wheel_zoom,hover,tap,reset,save"

    colors = get_bokeh_colors()
    df=df.rename(columns={'level_0':'predictor'})
    df['color'] = [colors[x] for x in df['predictor']]
    p = figure(title = "Summary", tools=TOOLS, width=500, height=500)
    p.xaxis.axis_label = 'binder_density'
    p.yaxis.axis_label = 'binders'

    #make metric for point sizes
    #df['point_size'] = df.binder_density
    source = ColumnDataSource(data=df)

    p.circle(x='binder_density', y='binders', line_color='black', fill_color='color',
             fill_alpha=0.4, size=10, source=source, legend_group='predictor')
    hover = p.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("name", "@name"),
        ("length", "@length"),
        ("binders", "@binders"),
        ("binder_density", "@binder_density"),
        ("top_peptide", "@top_peptide"),
        ("max_score", "@max_score"),
    ])
    p.toolbar.logo = None
    if savepath != None:
        url = "http://localhost:8000/sequence?savepath=%s&name=@name" %savepath
        taptool = p.select(type=TapTool)
        taptool.callback = OpenURL(url=url)

    callback = CustomJS(args=dict(source=source), code="""
        var data = source.data;
        var f = cb_obj.value

        data['x'] = f
        source.trigger('change');
        source.change.emit();
    """)

    from bokeh.layouts import widgetbox
    from bokeh.models.widgets import Select
    menu = [(i,i) for i in df.columns]
    select = Select(title='X', value='A', options=list(df.columns), width=8)
    select.js_on_change('value', callback)

    #layout = column(p, select, sizing_mode='scale_width')
    return p

def bokeh_plot_tracks(preds, title='', n=2, name=None, cutoff=.95, cutoff_method='default',
                width=None, height=None, x_range=None, tools=True, palette='Set1',
                seqdepot=None, exp=None):
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
    from bokeh.models import Range1d, HoverTool, FactorRange, ColumnDataSource, Text, Rect
    from bokeh.plotting import figure

    if tools == True:
        tools="xpan, xwheel_zoom, hover, reset, save"
    else:
        tools=''
    if width == None:
        width=1000
        sizing_mode='scale_width'
    else:
        sizing_mode='fixed'
    alls=1
    seqlen=0
    for P in preds:
        if P.data is None or len(P.data)==0:
            continue
        seqlen = get_seq_from_binders(P, name=name)
        #print (seqlen)
        alls += len(P.data.groupby('allele'))
    if seqlen == 0:
        return
    if height==None:
        height = 140+10*alls
    if x_range == None:
        x_range = Range1d(0, seqlen, bounds='auto')
    yrange = Range1d(start=0, end=alls+3)

    plot = figure(title=title, plot_width=width, sizing_mode=sizing_mode,
                    plot_height=height, y_range=yrange, x_range=x_range,
                    y_axis_label='allele',
                    tools=tools)
    h=3

    if exp is not None:
        plotExp(plot, exp)

    colors = get_bokeh_colors(palette)
    x=[];y=[];allele=[];widths=[];clrs=[];peptide=[]
    predictor=[];position=[];score=[];leg=[];seqs=[];text=[]
    l=80
    i=0
    for pred in preds:
        m = pred.name
        df = pred.data
        seq = base.sequence_from_peptides(df)
        if df is None or len(df) == 0:
            print('no data to plot for %s' %m)
            continue
        if name != None:
            df = df[df.name==name]

        sckey = pred.scorekey
        binders = pred.get_binders(name=name, cutoff=cutoff, cutoff_method=cutoff_method)
        #print (cutoff, n)
        pb = pred.promiscuous_binders(n=n, name=name, cutoff=cutoff, cutoff_method=cutoff_method)
        if len(pb) == 0:
            continue
        l = base.get_length(pb)
        grps = df.groupby('allele')
        alleles = grps.groups.keys()
        #seqs.extend([seq for i in alleles])
        #t = [i for s in list(seqs) for i in s]
        #text.extend(t)
        if len(pb)==0:
            continue
        c = colors[m]

        leg.append(m)
        seqlen = df.pos.max()+l

        for a,g in grps:
            b = binders[binders.allele==a]
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
        i+=1

    data = dict(x=x,y=y,allele=allele,peptide=peptide,width=widths,color=clrs,
                predictor=predictor,position=position,score=score)
    source = ColumnDataSource(data=data)
    plot.rect(x='x', y='y', source=source, width='width', height=0.8,
             legend_group='predictor',
             color='color',line_color='gray',alpha=0.7)

    #glyph = Text(x="x", y="y", text="text", text_align='center', text_color="black",
    #             text_font="monospace", text_font_size="10pt")
    #plot.add_glyph(source, glyph)

    hover = plot.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("allele", "@allele"),
        ("position", "@position"),
        ("peptide", "@peptide"),
        ("score", "@score"),
        ("predictor", "@predictor"),
    ])

    plot.xaxis.major_label_text_font_size = "9pt"
    plot.xaxis.major_label_text_font_style = "bold"
    plot.ygrid.grid_line_color = None
    plot.xgrid.minor_grid_line_alpha = 0.1
    plot.xgrid.minor_grid_line_color = 'gray'
    #plot.xgrid.minor_grid_line_dash = [6, 4]
    plot.yaxis.major_label_text_font_size = '0pt'
    #plot.xaxis.major_label_orientation = np.pi/4
    plot.min_border = 10
    plot.background_fill_color = "#fafaf4"
    plot.background_fill_alpha = 0.5
    plot.legend.orientation = "horizontal"
    plot.legend.location = "bottom_right"
    #plot.legend.label_text_font_size = 12
    plot.toolbar.logo = None
    plot.toolbar_location = "right"
    return plot

def bokeh_plot_sequence(preds, name=None, n=2, cutoff=.95, cutoff_method='default',
                        width=1000, color_sequence=False, title=''):
    """Plot sequence view of binders """

    from bokeh.plotting import figure
    from bokeh.models import ColumnDataSource, LinearAxis, Grid, Range1d, Text, Rect, CustomJS, Slider, RangeSlider, FactorRange
    from bokeh.layouts import gridplot, column

    colors = []
    seqs = []
    text = []
    alleles = []
    ylabels = []
    pcolors = get_bokeh_colors()

    for P in preds:
        print (P.name)
        df = P.data
        #get sequence from results dataframe
        seq = base.sequence_from_peptides(df)
        l = base.get_length(df)
        b = P.get_binders(name=name, cutoff=cutoff, cutoff_method=cutoff_method)
        pb = P.promiscuous_binders(name=name, cutoff=cutoff, n=n, cutoff_method=cutoff_method)
        b = b[b.pos.isin(pb.pos)] #only promiscuous

        grps = b.groupby('allele')
        al = list(grps.groups)
        alleles.extend(al)
        ylabels.extend([P.name+' '+i for i in al])
        currseq=[seq for i in al]
        seqs.extend(currseq)
        t = [i for s in currseq for i in s]
        text.extend(t)
        print (len(seqs),len(text))
        for a in al:
            pos=[]
            f = list(b[b.allele==a].pos)
            for i in f:
                pos.extend(np.arange(i,i+l))
            if color_sequence is True:
                c = plotting.get_sequence_colors(seq)
            else:
                c = ['white' for i in seq]

            for i in pos:
                c[i] = pcolors[P.name]
            colors.extend(c)

    #put into columndatasource for plotting
    N = len(seqs[0])
    S = len(alleles)
    x = np.arange(1, N+1)
    y = np.arange(0,S,1)
    xx, yy = np.meshgrid(x, y)
    gx = xx.ravel()
    gy = yy.flatten()
    recty = gy+.5

    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs)*15+60
    x_range = Range1d(0,N+1, bounds='auto')
    L=100
    if len(seq)<100:
        L=len(seq)
    view_range = (0,L)
    viewlen = view_range[1]-view_range[0]

    fontsize="8.5pt"
    tools="xpan, reset, save"
    p = figure(title=title, plot_width=width, plot_height=plot_height, x_range=view_range, y_range=ylabels, tools=tools,
               min_border=0, sizing_mode='stretch_both', lod_factor=10,  lod_threshold=1000)
    seqtext = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                 text_font="monospace", text_font_size=fontsize)
    rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors", line_color='gray', fill_alpha=0.6)

    p.add_glyph(source, rects)
    p.add_glyph(source, seqtext)
    p.xaxis.major_label_text_font_style = "bold"
    p.grid.visible = False
    p.toolbar.logo = None

    #preview view (no text)
    p1 = figure(title=None, plot_width=width, plot_height=S*3+5, x_range=x_range, y_range=(0,S), tools=[],
                    min_border=0, sizing_mode='stretch_width', lod_factor=10, lod_threshold=10)
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors", line_color=None, fill_alpha=0.6)
    previewrect = Rect(x=viewlen/2,y=S/2, width=viewlen, height=S*.99, line_color='darkblue', fill_color=None)
    p1.add_glyph(source, rects)
    p1.add_glyph(source, previewrect)
    p1.yaxis.visible = False
    p1.grid.visible = False
    p1.toolbar_location = None

    #callback for slider move
    jscode="""
        var start = cb_obj.value[0];
        var end = cb_obj.value[1];
        x_range.setv({"start": start, "end": end})
        rect.width = end-start;
        rect.x = start+rect.width/2;
        var fac = rect.width/width;
        console.log(fac);
        if (fac>=.14) { fontsize = 0;}
        else { fontsize = 8.5; }
        text.text_font_size=fontsize+"pt";
    """
    callback = CustomJS(
        args=dict(x_range=p.x_range,rect=previewrect,text=seqtext,width=p.plot_width), code=jscode)
    slider = RangeSlider (start=0, end=N, value=(0,L), step=10)#, callback_policy="throttle")
    slider.js_on_change('value_throttled', callback)

    #callback for plot drag
    jscode="""
        start = parseInt(range.start);
        end = parseInt(range.end);
        slider.value[0] = start;
        rect.width = end-start;
        rect.x = start+rect.width/2;
    """
    p.x_range.callback = CustomJS(args=dict(slider=slider, range=p.x_range, rect=previewrect),
                                  code=jscode)

    p = gridplot([[p1],[p],[slider]], toolbar_location="below", merge_tools=False)
    return p

def bokeh_plot_grid(pred, name=None, width=None, palette='Blues', **kwargs):
    """Plot heatmap of binding results for a predictor."""

    from bokeh.plotting import figure
    from bokeh.models import (Range1d,HoverTool,FactorRange,ColumnDataSource,
                              LinearColorMapper,LogColorMapper,callbacks,DataRange)
    from bokeh.palettes import all_palettes
    TOOLS = "xpan, xwheel_zoom, hover, reset, save"

    if width == None:
        sizing_mode = 'scale_width'
        width=900
    else:
        sizing_mode = 'fixed'
    P=pred
    df = P.data
    if df is None:
        return

    cols = ['allele','pos','peptide',P.scorekey]
    #d = df[cols].copy()
    b = P.get_binders(name=name,**kwargs)
    d = P.data.copy()
    #mark binders
    mask = d.index.isin(b.index)
    d['binder'] = mask

    l = base.get_length(df)
    grps = df.groupby('allele')
    alleles = grps.groups
    seqlen = get_seq_from_binders(P, name)
    seq = base.seq_from_binders(df)
    height = 300
    alls = len(alleles)
    x_range = Range1d(0,seqlen-l+1, bounds='auto')
    #x_range = list(seq)
    y_range = df.allele.unique()
    val = P.scorekey
    cut = P.cutoff
    if P.name not in ['tepitope']:
        d['score1'] = d[val].apply( lambda x: 1-math.log(x, 50000))
        val='score1'
    d[val][d.binder==False] = min(d[val])

    source = ColumnDataSource(d)
    colors = all_palettes[palette][7]
    mapper = LinearColorMapper(palette=colors, low=d[val].max(), high=d[val].min())

    p = figure(title=P.name+' '+name,
               x_range=x_range, y_range=y_range,
               x_axis_location="above", plot_width=width, plot_height=height,
               tools=TOOLS, toolbar_location='below', sizing_mode=sizing_mode)
    p.rect(x="pos", y="allele", width=1, height=1,
           source=source,
           fill_color={'field': val,'transform':mapper},
           line_color='gray', line_width=.1)
    p.select_one(HoverTool).tooltips = [
         ('allele', '@allele'),
         (P.scorekey, '@%s{1.11}' %P.scorekey),
         ('pos', '@pos'),
         ('peptide', '@peptide')
    ]
    p.toolbar.logo = None
    p.yaxis.major_label_text_font_size = "10pt"
    p.yaxis.major_label_text_font_style = "bold"
    return p

def bokeh_plot_bar(preds, name=None, allele=None, title='', width=None, height=100,
                    palette='Set1', tools=True, x_range=None):
    """Plot bars combining one or more prediction results for a set of
    peptides in a protein/sequence"""

    from bokeh.models import Range1d,HoverTool,ColumnDataSource
    from bokeh.plotting import figure
    from bokeh.transform import dodge
    from bokeh.core.properties import value

    height = 180
    seqlen = 0
    if width == None:
        width=700
        sizing_mode='scale_width'
    for P in preds:
        if P.data is None or len(P.data)==0:
            continue
        seqlen = get_seq_from_binders(P)

    if x_range == None:
        x_range = Range1d(0,seqlen)
    y_range = Range1d(start=0, end=1)
    if tools == True:
        tools="xpan, xwheel_zoom, reset, hover"
    else:
        tools=None
    plot = figure(title=title,plot_width=width,sizing_mode=sizing_mode,
                    plot_height=height, y_range=y_range, x_range=x_range,
                    y_axis_label='rank',
                    tools=tools)
    colors = get_bokeh_colors(palette)
    data = {}
    mlist = []
    for pred in preds:
        m = pred.name
        df = pred.data
        if df is None or len(df) == 0:
            continue
        if name != None:
            df = df[df.name==name]
        grps = df.groupby('allele')
        alleles = grps.groups.keys()
        if allele not in alleles:
            continue
        #print (m, alleles, allele)
        df = df[df.allele==allele]
        df = df.sort_values('pos').set_index('pos')
        key = pred.scorekey
        X = df[key]
        X = (X+abs(X.min())) / (X.max() - X.min())
        data[m] = X.values
        data['pos'] = list(X.index)
        #data['peptide'] = df.peptide.values
        mlist.append(m)
    source = ColumnDataSource(data)

    w = round(1.0/len(mlist),1)-.1
    i=-w/2
    for m in mlist:
        #m = pred.name
        c = colors[m]
        plot.vbar(x=dodge('pos', i, range=plot.x_range), top=m, width=w, source=source,
                   color=c, legend=value(m), alpha=.8)
        i+=w

    hover = plot.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        #("allele", "@allele"),
        ("pos", "@pos") ])

    plot.min_border = 10
    plot.background_fill_color = "beige"
    plot.background_fill_alpha = 0.5
    plot.toolbar.logo = None
    plot.toolbar_location = "right"
    plot.legend.location = "top_right"
    plot.legend.orientation = "horizontal"
    return plot

def bokeh_vbar(x, height=200, title='', color='navy'):

    from bokeh.plotting import figure
    from bokeh.models import ColumnDataSource

    source = ColumnDataSource(data={'chr':list(x.index),'x':range(len(x)),'y':x.values})
    plot = figure(title=title, x_range = list(x.index), plot_height=height, tools='save,reset')
    plot.vbar(x='chr',top='y', width=.8, bottom=0,source=source, color=color)
    plot.ygrid.grid_line_color = None
    plot.xgrid.grid_line_color = None
    plot.xaxis.major_label_orientation = np.pi/4
    return plot

def bokeh_pie_chart(df, title='', radius=.5, width=400, height=400, palette='Spectral'):
    """Bokeh pie chart"""

    from bokeh.plotting import figure
    from bokeh.models import HoverTool,ColumnDataSource
    from math import pi

    s = df.cumsum()/df.sum()
    cats = s.index
    p=[0]+list(s)
    #print (p)
    starts = [1/2*pi-(i*2*pi) for i in p[:-1]]
    ends = [1/2*pi-(i*2*pi) for i in p[1:]]

    from bokeh.palettes import brewer
    n = len(s)
    pal = brewer[palette][6]
    source = ColumnDataSource(
                dict(x=[0 for x in s], y=[0 for x in s],
                  radius = [radius for x in s],
                  category= cats,
                  starts=starts,
                  ends=ends,
                  colors=pal,
                  counts = df
                  ))
    plot = figure(title=title, plot_width=width, plot_height=height, tools='save,reset')
    plot.wedge(x='x', y='y', radius='radius', direction="clock", fill_color='colors', color='black',
                start_angle='starts', end_angle='ends', legend='category', source=source)
    plot.axis.visible = False
    plot.ygrid.visible = False
    plot.xgrid.visible = False
    #hover = plot.select(dict(type=HoverTool))
    #hover.tooltips = [
    #    ('category', '@category'),
    #    ('percents','@counts')
    #]
    return plot

def plot_tracks(preds, name, n=1, cutoff=.95, cutoff_method='default', regions=None,
                legend=False, colormap='Paired', figsize=None, ax=None, **kwargs):
    """
    Plot binders as bars per allele using matplotlib.
    Args:
        preds: list of one or more predictors
        name: name of protein to plot
        n: number of alleles binder should be found in to be displayed
        cutoff: percentile cutoff to determine binders to show

    """

    import matplotlib as mpl
    import pylab as plt
    from matplotlib.patches import Rectangle
    if ax == None:
        if figsize==None:
            h = sum([len(p.data.groupby('allele')) for p in preds])
            w = 10
            h = round(h*.1+2)
            figsize = (w,h)
        #plt.clf()
        fig = plt.figure(figsize=figsize,facecolor='white')
        ax = fig.add_subplot(111)

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
            if name not in df.name.unique():
                print ('no such sequence %s' %name)
                continue
            df = df[df.name==name]
        sckey = pred.scorekey

        binders = pred.get_binders(name=name, cutoff=cutoff,
                                   cutoff_method=cutoff_method)
        #print (binders)
        pb = pred.promiscuous_binders(binders=binders, n=n)

        if len(pb) == 0:
            continue
        l = base.get_length(pb)
        seqlen = df.pos.max()+l
        #print (name,m,df.pos.max(),l,seqlen)
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
                rect = ax.add_patch(Rectangle((x,y), l, 1, facecolor=c, edgecolor='black',
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
    if regions is not None:
        r = regions[regions.name==name]
        coords = (list(r.start),list(r.end-r.start))
        coords = zip(*coords)
        plot_regions(coords, ax, color='gray')
    if legend == True:
        ax.legend(handles, leg, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=3)
    plt.tight_layout()
    return ax

def plot_regions(coords, ax, color='red', label='', alpha=0.6):
    """Highlight regions in a prot binder plot"""

    from matplotlib.patches import Rectangle
    #l = len(seqs.head(1)['key'].max())
    h = ax.get_ylim()[1]
    for c in coords:
        x,l = c
        ax.add_patch(Rectangle((x,0), l, h,
            facecolor=color, lw=.8, alpha=alpha, zorder=0))
    return

def draw_labels(labels, coords, ax):
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

def plot_bars(P, name, chunks=1, how='median', cutoff=20, color='black'):
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
    l= base.get_length(df)
    seqlen = df.pos.max()+l

    funcs = {'median': np.median, 'mean': np.mean, 'sum': np.sum}
    grps = df.groupby('pos')
    key = P.scorekey
    X = grps.agg({key: np.median, 'peptide': base.first})
    q = (1-cutoff/100.) #score quantile value
    cutoff = X[key].quantile(q)
    X[key][X[key]<cutoff] = np.nan
    if len(X)<20:
        chunks = 1
    seqlist = X.peptide.apply( lambda x : x[0])
    seqchunks = np.array_split(X.index, chunks)

    f,axs = plt.subplots(chunks,1,figsize=(15,2+2.5*chunks))
    if chunks == 1:
        axs = [axs]
    else:
        axs = list(axs.flat)
    for c in range(chunks):
        #print (c)
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

def plot_bcell(plot,pred,height,ax=None):
    """Line plot of iedb bcell results"""

    x = pred.data.Position
    y = pred.data.Score
    h = height
    y = y+abs(min(y))
    y = y*(h/max(y))+3
    #plot.line(x, y, line_color="red", line_width=2, alpha=0.6,legend='bcell')
    ax.plot(x,y,color='blue')
    return

def plot_seqdepot(annotation, ax):
    """Plot sedepot annotations - replace with generic plot coords track"""

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
        #print (pos,widths,text)
        bbox_args = dict(boxstyle='round', fc='white', lw=1, alpha=0.8)
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
            plot_regions(coords, ax, color='gray')
        #labels = list(r.peptide)
        #plotting.draw_labels(labels, coords, ax)
        if genome is not None:
            p = genome[genome['locus_tag']==prot]
            seq = p.translation.iloc[0]
            from . import analysis
            sd = analysis.get_seqdepot(seq)['t']
            plot_seqdepot(sd, ax)
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

    import pylab as plt
    import seaborn as sns
    df = P.data[P.data.name==name].sort_values('pos')
    w = 10
    l= base.get_length(df)
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
        vmin=min(X.min()); vmax=max(X.max())
        center = vmin+(vmax-vmin)/2
        sns.heatmap(p, cmap=cmap, cbar_kws={"shrink": .5},
                    vmin=vmin, vmax=vmax, #center=center,
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

def binders_to_coords(df):
    """Convert binder results to dict of coords for plotting"""

    coords = {}
    if not 'start' in df.columns:
        df=base.get_coords(df)
    if 'start' in df.columns:
        for i,g in df.groupby('name'):
            l = g.end-g.start
            coords[i] = zip(g.start,l)
    return coords

def plot_overview(genome, coords=None, cols=2, colormap='Paired',
                  legend=True, figsize=None):
    """
    Plot regions of interest in a group of protein sequences. Useful for
    seeing how your binders/epitopes are distributed in a small genome or subset of genes.
    Args:
        genome: dataframe with protein sequences
        coords: a list/dict of tuple lists of the form {protein name: [(start,length)..]}
        cols: number of columns for plot, integer
    """

    import pylab as plt
    if type(coords) is list:
        coords = { i:coords[i] for i in range(len(coords)) }
        legend=False
    import matplotlib as mpl
    import seaborn as sns
    #sns.reset_orig()
    cmap = mpl.cm.get_cmap(colormap)
    t = len(coords)
    colors = [cmap(float(i)/t) for i in range(t)]
    from matplotlib.patches import Rectangle

    names = [coords[c].keys() for c in coords][0]

    df = genome[genome.locus_tag.isin(names)]
    rows = int(np.ceil(len(names)/float(cols)))
    if figsize == None:
        h = round(len(names)*.2+10./cols)
        figsize = (14,h)
    f,axs=plt.subplots(rows,cols,figsize=figsize)
    grid=axs.flat
    rects = {}
    i=0
    for idx,prot in df.iterrows():
        ax=grid[i]
        protname = prot.locus_tag
        seq = prot.translation
        if 'description' in prot:
            title = prot.description
        else:
            title = protname
        y=0
        for label in coords:
            c = coords[label]
            if not protname in c:
                continue
            vals = c[protname]
            #print vals
            for v in vals:
                x,l = v[0],v[1]
                rect = ax.add_patch(Rectangle((x,y), l, .9,
                                facecolor=colors[y],
                                lw=1.2, alpha=0.8))
                if len(v)>2:
                    s = v[2]
                    bbox_args = dict(fc=colors[y], lw=1.2, alpha=0.8)
                    ax.annotate(s, (x+l/2, y),
                        fontsize=12, ha='center', va='bottom')
            if not label in rects:
                rects[label] = rect
            y+=1
        i+=1
        slen = len(seq)
        w = round(float(slen)/20.)
        w = math.ceil(w/20)*20
        ax.set_xlim(0, slen)
        ax.set_ylim(0, t)
        ax.set_xticks(np.arange(0, slen, w))
        ax.set_yticks([])
        ax.set_title(title, fontsize=16, loc='right')

    if i|2!=0 and cols>1:
        try:
            f.delaxes(grid[i])
        except:
            pass
    if legend == True:
        f.legend(rects.values(), rects.keys(), loc=4)
    plt.tight_layout()
    return

def seqdepot_to_coords(sd, key='pfam27'):
    """
    Convert seqdepot annotations to coords for plotting
    """

    coords=[]
    if len(sd['t'])==0 or not key in sd['t']:
        return []
    x = sd['t'][key]
    #print x

    if key in ['pfam27','pfam28']:
        coords = [(i[1],i[2]-i[1],i[0]) for i in x]
    elif key in ['gene3d','prints']:
        coords = [(i[2],i[3]-i[2],i[1]) for i in x]
    elif key == 'tmhmm':
        coords = [(i[0],i[1]-i[0]) for i in x]
    elif key == 'signalp':
        x = x.items()
        coords = [(i[1],10,'SP') for i in x]
    return coords

def get_seqdepot_annotation(genome, key='pfam27'):
    """
    Get seqdepot annotations for a set of proteins in dataframe.
    """
    from . import seqdepot
    annot={}
    for i,row in genome.iterrows():
        n = row.locus_tag
        seq = row.translation
        #print n,seq
        sd = seqdepot.new()
        aseqid = sd.aseqIdFromSequence(seq)
        result = sd.findOne(aseqid)
        #for x in result['t']:
        #    print x, result['t'][x]
        x = seqdepot_to_coords(result, key)
        annot[n] = x
    return annot
