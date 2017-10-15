#!/usr/bin/env python

"""
    epitopepredict, methods for supporting web app
    Created Sep 2017
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys,os,glob
import pandas as pd
import numpy as np
from . import base, plotting, sequtils
from bokeh.models import ColumnDataSource, Slider
from bokeh.models.widgets import DataTable, TableColumn, Select, Button, Slider, TextInput
from bokeh.layouts import row, column, gridplot, widgetbox, layout
from bokeh.embed import components

path = 'results'
predictors = base.predictors
wikipage = 'https://github.com/dmnfarrell/epitopepredict/wiki'
plotkinds = ['tracks','bar','text']

def get_file_lists(path):

    names = []
    for p in predictors:
        files = glob.glob(os.path.join(path, p, '*.csv'))
        n = [os.path.splitext(os.path.basename(i))[0] for i in files]
        names.extend(n)
    names = set(names)
    names = sorted(names)
    return names

def get_results(path, predictor, name):

    filename = os.path.join(path, predictor, name)
    P = base.get_predictor(predictor)
    P.load(filename+'.csv')
    #print filename
    #print P.data
    return P

def get_results_info(P):
    """Info on sequence used for prediction"""

    df = P.data
    if df is None:
        return ''
    df = df.drop_duplicates('pos').sort_values('pos')
    #l = base.get_length(df)
    seq = sequence_from_peptides(df)
    l = len(seq)
    return {'length':l}

def sequence_from_peptides(df):
    x = df.peptide.str[0]
    x = ''.join(x)
    return x

def get_predictors(path, name):
    """Get a set of predictors with available results"""

    if name==None:
        return []
    preds = []
    for pred in predictors:
        P = get_results(path, pred, name)
        preds.append(P)
    return preds

def create_figures(preds, name='', kind='tracks', cutoff=5, n=2,
                   cutoff_method='default', **kwargs):
    """Get plots of binders for single protein/sequence"""

    figures = []
    if kind == 'tracks':
        plot = plotting.bokeh_plot_tracks(preds, title=name,
                         width=700, palette='Set1', cutoff=float(cutoff), n=int(n),
                         cutoff_method=cutoff_method)

    elif kind == 'bar':
        plot = plotting.bokeh_plot_bar(preds, title=name, width=700 )
    if plot is not None:
        figures.append(plot)
    return figures

def create_bokeh_table(path, name):
    """Create table of prediction data"""

    P = get_results(path, 'tepitope', name)
    if P.data is None:
        return
    df = P.data[:10]
    data = dict(
        peptide=df.peptide.values,
        pos=df.pos.values,
        score=df.score.values,
        allele=df.allele.values
    )
    #print (df)
    source = ColumnDataSource(data)
    columns = [
            TableColumn(field="peptide", title="peptide"),
            TableColumn(field="pos", title="pos"),
            TableColumn(field="score", title="score"),
            TableColumn(field="allele", title="allele"),
        ]
    table = DataTable(source=source, columns=columns, width=400, height=280)
    return table

def create_binder_tables(preds, name, promiscuous=False,
                         escape=False, classes='', **kwargs):
    """Create table of prediction data"""

    tables = {}
    import pylab as plt
    cm = plt.get_cmap('Reds')
    for P in preds:
        df = P.getBinders(name=name, **kwargs)
        if df is None:
            continue
        df = df.reset_index(drop=True)
        s = df.style\
              .set_table_attributes('class="%s"' %classes)
              #.background_gradient(subset=[P.scorekey], cmap=cm) #%classes
        #t = df.to_html(escape=escape, classes=classes, index=True)
        t = s.render(index=False)
        tables[P.name] = t
    return tables

def dict_to_html(data):

    s = ''
    for k in data:
        s += '<a>%s: %s</a><br>' %(k,data[k])

    return s

def tabbed_html(items):
    """Create html for a set of tabbed divs from dict of html code, one for
       each tab. Uses css classes defined in static/custom.css"""

    name = 'tab-group'
    html = '<div class="tabs">\n'
    for t in items:
        html += '<div class="tab">\n'
        html += '<input type="radio" id="%s" name="%s" checked>\n' %(t,name)
        html +=	'<label for="%s">%s</label>\n' %(t,t)
        html +=	'<div class="content">\n'
        html += items[t]
        #html += '<p>%s</p>' %t
        html +=	'</div></div>\n'
    html += '</div>'
    print (html)
    return html

def create_widgets():

    select = Select(title="Name:", value="name", options=["foo", "bar"])
    slider = Slider(start=0, end=100, value=5, step=.5, title="Cutoff")
    button = Button(label="Submit", button_type="success")
    #text_input = TextInput(value="default", title="Name:")
    return widgetbox([button,select,slider], width=200)

def test():
    from bokeh.io import output_file, show
    path = 'results'
    name = 'Rv0011c'
    kwargs ={'cutoff_method':'default'}
    preds = get_predictors(path, name)
    plots = create_figures(preds)
    #table = create_bokeh_table(path, name)
    tables = create_binder_tables(preds, name)
    grid = gridplot(plots, ncols=1, merge_tools=True)
    widgets = create_widgets()
    l = layout([[ plots, widgets ]], ncols=2, nrows=1)
    #script, div = components(l)
    show(l)

if __name__ == "__main__":
    test()