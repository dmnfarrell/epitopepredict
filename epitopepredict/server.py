#!/usr/bin/env python

"""
    epitopepredict flask server app for viewing results
    Created August 2017
    Copyright (C) Damien Farrell
"""

import os,sys,glob
from flask import Flask, render_template, request
from wtforms import Form, TextField, validators, StringField, SelectField, FloatField
from wtforms import FileField, SubmitField
import pandas as pd
import numpy as np
#from bokeh.Charts import Histogram, Bar, Scatter
from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.layouts import row, column, gridplot
from bokeh.models.widgets import Button, RadioButtonGroup, Select, Slider

from epitopepredict import base, plotting

app = Flask(__name__)
path = 'results'
predictors = base.predictors

class ControlsForm(Form):
    name = SelectField('name', choices=[])
    path = TextField('path', default='results')
    cutoff = FloatField('cutoff', default=5)
    n = TextField('n', default='2')
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

def get_seq_info(P, name):
    l = base.get_length(P.data)
    return {'n-mer length':l}

# Create the main plot

def create_figure(name=None, kind='tracks', cutoff=5, n=2):
    """Get plot of binders for single protein/sequence"""

    if name==None:
        return []
    figures = []
    preds = []
    for pred in predictors:
        P = get_results(path, pred, name)
        preds.append(P)
    plot = plotting.bokeh_plot_tracks(preds, title=name,
                        width=800, palette='Set1', cutoff=float(cutoff), n=int(n))

    if plot is not None:
        figures.append(plot)
    return figures

def create_figures(name=None, kind='tracks'):
    """Create multiple separate figures"""

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
        plot = plotting.bokeh_plot_tracks([P], title=pred+' '+name, x_range=xr,
                            width=800, height=180, palette='Set1')

        if plot is not None:
            figures.append(plot)
    return figures

@app.route('/')
def index():
    """index page"""

    path = request.args.get("path")
    if path == None: path= 'results'
    names = get_file_lists(path)
    current_name = request.args.get("name")
    if current_name is None: current_name='Rv0011c'
    cutoff = request.args.get("cutoff")
    if cutoff is None: cutoff=5
    n = request.args.get("n")
    if n is None: n=2
    print cutoff

    form = ControlsForm()
    form.path.data = path
    form.name.choices = [(i,i) for i in names]
    form.name.data = current_name
    form.cutoff.data = cutoff
    form.n.data = n

    #get_seq_info(P, name)
    plots = create_figure(current_name, cutoff=cutoff, n=n)
    if len(plots) > 0:
        grid = gridplot(plots, ncols=1, merge_tools=True, #sizing_mode='stretch_both',
                        toolbar_options=dict(logo='grey'))
        script, div = components(grid)
    else:
        script=''; div=''
    return render_template("index.html", form=form, script=script, div=div,
            path=path, names=names, current_name=current_name)

if __name__ == '__main__':
	app.run(port=5000, debug=True)
