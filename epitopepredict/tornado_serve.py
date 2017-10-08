#!/usr/bin/env python

"""
    epitopepredict server app for viewing results, uses tornado
    Created Sep 2017
    Copyright (C) Damien Farrell
"""

import sys,os,glob
import pandas as pd
import numpy as np
from epitopepredict import base, plotting

import tornado.ioloop
import tornado.web
from wtforms_tornado import Form
from wtforms import TextField, StringField, SelectField, FloatField
from tornado.web import RequestHandler
from bokeh.util.browser import view
from bokeh.models import ColumnDataSource, Slider
from bokeh.plotting import figure
from bokeh.layouts import row, column, gridplot, widgetbox
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

def get_seq_info(P):
    df = P.data
    l = base.get_length(df)
    seq = sequence_from_peptides(df)
    return {'n-mer':l, 'sequence':seq}

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

class ControlsForm(Form):
    name = SelectField('name', choices=[])
    path = TextField('path', default='results')
    cutoff = FloatField('cutoff', default=5)
    n = TextField('n', default='2')
    kinds = [(i,i) for i in plotkinds]
    kind = SelectField('plot kind', choices=kinds)
    #submit = SubmitField()

class MainHandler(RequestHandler):
    """Handler for main results page"""
    def get(self):

        args = self.request.arguments
        form = ControlsForm()
        print args
        if 'path' in args:
           path = args['path'][0]
        else:
            path = 'results'
        if not os.path.exists(path):
            msg = help_msg()
            self.render('index.html', script='', div='', form=form, msg=msg)
            return
        if 'name' in args:
            current_name = args['name'][0]
        else:
            current_name = 'Rv0011c'
        cutoff=5
        n=2
        kind='tracks'
        names = get_file_lists(path)

        form.path.data = path
        form.name.choices = [(i,i) for i in names]
        form.name.data = current_name

        preds = get_predictors(path, current_name)
        plots = create_figures(preds, current_name, cutoff=cutoff, n=n, kind=kind)
        #info = get_seq_info(preds[0])['sequence']

        if len(plots) > 0:
            grid = gridplot(plots, ncols=1, merge_tools=True, #sizing_mode='stretch_both',
                            toolbar_options=dict(logo=None))
        script, div = components(grid)

        text = '<h4> using tornado.... </h4>'
        self.render('index.html', script=script, div=div, form=form, msg='')

    @staticmethod
    def modify_doc(doc):
        return

settings = dict(
        template_path=os.path.join(os.path.dirname(__file__), "templates"),
        static_path=os.path.join(os.path.dirname(__file__), "static"),
        autoescape=None,
        xsrf_cookies=True,
        debug=True)

def main():
    #bokeh_app = Application(FunctionHandler(IndexHandler.modify_doc))
    #bokeh_server = Server({'/main': bokeh_app},
    #                      io_loop=io_loop,
    #                      extra_patterns=[('/', IndexHandler)],
    #                      allow_websocket_origin=['localhost:5006'])
    #bokeh_server.start()
    handlers = [ (r"/", MainHandler) ]
    app = tornado.web.Application(handlers, **settings)
    #app.listen(8888)
    http_server = tornado.httpserver.HTTPServer(app)
    http_server.listen(8888)
    io_loop = tornado.ioloop.IOLoop.current()
    #io_loop.add_callback(view, "http://localhost:8888/")
    view("http://localhost:8888/")
    io_loop.start()


if __name__ == "__main__":
    main()