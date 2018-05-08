#!/usr/bin/env python

"""
    MHCpredict plotting
    Created February 2016
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, math
from collections import OrderedDict
import matplotlib
matplotlib.use('agg', warn=False)
import pylab as plt
import numpy as np
import pandas as pd
from . import base

colormaps={'tepitope':'Greens','netmhciipan':'Oranges','iedbmhc2':'Pinks',
               'threading':'Purples','iedbmhc1':'Blues'}
defaultcolors = {'tepitope':'green','netmhciipan':'orange',
           'iedbmhc1':'blue','iedbmhc2':'pink','threading':'purple'}

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
    if P.data is None or len(P.data)==0:
        return
    l = len(P.data.iloc[0].peptide)
    seqlen = P.data.pos.max()+l
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

def bokeh_test(height=400):

    from bokeh.models import ColumnDataSource
    from bokeh.plotting import figure
    data = {'x_values': np.random.random(20),
            'y_values': np.random.random(20)}
    source = ColumnDataSource(data=data)
    p = figure(plot_height=height)
    p.circle(x='x_values', y='y_values', radius=.04, line_color='black', fill_alpha=.6, source=source)
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
             fill_alpha=0.4, size=10, source=source, legend='predictor')
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

def bokeh_plot_tracks(preds, title='', n=2, name=None, cutoff=5, cutoff_method='default',
                width=None, height=None, x_range=None, tools=True,
                palette='Set1',
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
    from bokeh.models import Range1d,HoverTool,FactorRange,ColumnDataSource
    from bokeh.plotting import figure

    if tools == True:
        tools="xpan, xwheel_zoom, hover, reset, save"
    else:
        tools=''
    if width == None:
        width=1000
        sizing_mode='scale_width'
    else:
        sizing_mode=None
    alls=1
    seqlen=0
    for P in preds:
        if P.data is None or len(P.data)==0:
            continue
        seqlen = get_seq_from_binders(P)
        alls += len(P.data.groupby('allele'))
    if seqlen == 0:
        return
    if height==None:
        height = 140+10*alls
    if x_range == None:
        x_range = Range1d(0,seqlen)
    yrange = Range1d(start=0, end=alls+3)

    plot = figure(title=title, plot_width=width, sizing_mode=sizing_mode,
                    plot_height=height, y_range=yrange, x_range=x_range,
                    y_axis_label='allele',
                    tools=tools)
    h=3
    #if bcell != None:
    #    plotBCell(plot, bcell, alls)
    if seqdepot != None:
        plotAnnotations(plot,seqdepot)
    if exp is not None:
        plotExp(plot, exp)

    colors = get_bokeh_colors(palette)
    x=[];y=[];allele=[];widths=[];clrs=[];peptide=[]
    predictor=[];position=[];score=[];leg=[]
    l=80
    i=0
    for pred in preds:
        m = pred.name
        df = pred.data
        if df is None or len(df) == 0:
            print('no data to plot for %s' %m)
            continue
        if name != None:
            df = df[df.name==name]
        sckey = pred.scorekey

        binders = pred.get_binders(name=name, cutoff=cutoff, cutoff_method=cutoff_method)
        #print (cutoff, n)
        pb = pred.promiscuous_binders(n=n, cutoff=cutoff, cutoff_method=cutoff_method)
        if len(pb) == 0:
            continue
        l = base.get_length(pb)
        grps = df.groupby('allele')
        alleles = grps.groups.keys()
        if len(pb)==0:
            continue
        c = colors[m]

        leg.append(m)
        seqlen = df.pos.max()+l

        for a,g in grps:
            #b = pred.getBinders(data=g, cutoff=cutoff)
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
    plot.rect(x='x',y='y', source=source, width='width', height=0.8,
             legend='predictor',
             color='color',line_color='gray',alpha=0.7)

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
    plot.yaxis.major_label_text_font_size = '0pt'
    #plot.xaxis.major_label_orientation = np.pi/4
    plot.min_border = 10
    plot.background_fill_color = "beige"
    plot.background_fill_alpha = 0.5
    plot.legend.orientation = "horizontal"
    plot.legend.location = "bottom_right"
    plot.toolbar.logo = None
    plot.toolbar_location = "right"
    return plot

def bokeh_plot_grid(pred, name, width=None, palette='Blues', **kwargs):
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
    b = P.get_binders(**kwargs)
    d = P.data.copy()
    #mark binders
    mask = d.index.isin(b.index)
    d['binder'] = mask

    l = base.get_length(df)
    grps = df.groupby('allele')
    alleles = grps.groups
    seqlen = get_seq_from_binders(P)
    seq = base.seq_from_binders(df)
    height = 300
    alls = len(alleles)
    x_range = Range1d(0,seqlen-l+1, bounds='auto')
    #x_range = list(seq)
    y_range = df.allele.unique()
    val=P.scorekey
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
           line_color='gray', line_width=.5)
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

def plot_tracks(preds, name, n=1, cutoff=5, value='score',
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
    if ax == None:
        if figsize==None:
            h = sum([len(p.data.groupby('allele')) for p in preds])
            w = 10
            h = round(h*.1+2)
            figsize = (w,h)
        plt.clf()
        fig = plt.figure(figsize=figsize)
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
            df = df[df.name==name]
        sckey = pred.scorekey

        binders = pred.get_binders(name=name, cutoff=cutoff, value=value)
        #binders = pred.binders
        #print (binders)
        pb = pred.promiscuous_binders(binders=binders, n=n, value=value)

        if len(pb) == 0:
            continue
        l = base.get_length(pb)
        seqlen = df.pos.max()+l
        print (name,m,df.pos.max(),l,seqlen)
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
    if legend == True:
        ax.legend(handles, leg, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=3)
    #plt.tight_layout()
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

    import seaborn as sns
    df = P.data[P.data.name==name].sort('pos')
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

    if type(coords) is list:
        coords = { i:coords[i] for i in range(len(coords)) }
        legend=False
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
        f.delaxes(grid[i])
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
