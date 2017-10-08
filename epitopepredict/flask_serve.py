#!/usr/bin/env python

"""
    epitopepredict flask server app for viewing results
    Created August 2017
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import os,sys,glob
from flask import Flask, render_template, request
from wtforms import Form, TextField, validators, StringField, SelectField, FloatField
from wtforms import FileField, SubmitField
import pandas as pd
import numpy as np

from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.layouts import row, column, gridplot, widgetbox
from bokeh.models.widgets import Button, RadioButtonGroup, Select, Slider
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn
from bokeh.models import ColumnDataSource

from epitopepredict import base, plotting

#apppath = os.path.dirname(os.path.abspath(__file__))
webapp = Flask(__name__)
path = 'results'
predictors = base.predictors
wikipage = 'https://github.com/dmnfarrell/epitopepredict/wiki'
plotkinds = ['tracks','bar','text']

class ControlsForm(Form):
    name = SelectField('name', choices=[])
    path = TextField('path', default='results')
    cutoff = FloatField('cutoff', default=5)
    n = TextField('n', default='2')
    kinds = [(i,i) for i in plotkinds]
    kind = SelectField('plot kind', choices=kinds)
    #submit = SubmitField()

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

def get_seq_info(P):
    df = P.data
    l = base.get_length(df)
    seq = sequence_from_peptides(df)
    return {'n-mer':l, 'sequence':seq}

def sequence_from_peptides(df):
    x = df.peptide.str[0]
    x = ''.join(x)
    return x

# Create the main plot

def get_predictors(name=None):
    """Get a set of predictors with available results"""

    if name==None:
        return []
    preds = []
    for pred in predictors:
        P = get_results(path, pred, name)
        preds.append(P)
    return preds

def create_figures(preds, name='', kind='tracks', cutoff=5, n=2):
    """Get plots of binders for single protein/sequence"""

    figures = []
    if kind == 'tracks':
        plot = plotting.bokeh_plot_tracks(preds, title=name,
                         width=800, palette='Set1', cutoff=float(cutoff), n=int(n))
    elif kind == 'bar':
        plot = plotting.bokeh_plot_bar(preds, title=name)
    if plot is not None:
        figures.append(plot)
    return figures

def create_multiple_figures(name=None, kind='bar'):
    """Create multiple separate figures for each predictor"""

    if name==None:
        return []
    figures = []
    plot = None
    for pred in predictors:
        P = get_results(path, pred, name)
        if plot is not None:
            xr = plot.x_range
        else:
            xr=None
        plot = plotting.bokeh_plot_bars([P], title=pred+' '+name, x_range=xr,
                            width=800, height=180)

        if plot is not None:
            figures.append(plot)
    return figures

def create_pred_table(path, name):
    """Create table of prediction data"""

    P = get_results(path, 'tepitope', name)
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

def help_msg():
    msg = '<a>results path not found, enter a folder with your results</a><br>'
    msg += '<a href="%s"> see help page</a>' %wikipage
    return msg

@webapp.route('/')
def index():
    """main index page"""

    path = request.args.get("path")
    if path == None: path= 'results'
    if not os.path.exists(path):
        msg = help_msg()
        return render_template("index.html",form=ControlsForm(),div='',msg=msg)

    names = get_file_lists(path)
    current_name = request.args.get("name")
    if current_name is None: current_name='Rv0011c'
    cutoff = request.args.get("cutoff")
    if cutoff is None: cutoff=5
    n = request.args.get("n")
    if n is None: n=2
    kind = request.args.get("kind")
    if kind is None: kind='tracks'

    form = ControlsForm()
    form.path.data = path
    form.name.choices = [(i,i) for i in names]
    form.name.data = current_name
    form.cutoff.data = cutoff
    form.n.data = n
    form.kind.data = kind

    script=''; div=''
    preds = get_predictors(current_name)
    plots = create_figures(preds, current_name, cutoff=cutoff, n=n, kind=kind)
    info = get_seq_info(preds[0])['sequence']

    if len(plots) > 0:
        grid = gridplot(plots, ncols=1, merge_tools=True, #sizing_mode='stretch_both',
                        toolbar_options=dict(logo=None))
        #script, div = components(grid)

    #table = widgetbox(create_pred_table(path, current_name))
    script, div = components({"plots": grid})#, "table": table})

    return render_template("index.html", form=form, script=script, div=div, msg='',
                           info=info, path=path, names=names, current_name=current_name)

def main():
    webapp.run(port=5000, debug=True)

if __name__ == '__main__':
	main()
