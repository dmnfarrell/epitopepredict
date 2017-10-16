#!/usr/bin/env python

"""
    epitopepredict server app for viewing results, uses tornado
    Created Sep 2017
    Copyright (C) Damien Farrell
"""

import sys,os,glob
import pandas as pd
import numpy as np
from epitopepredict import base, web

import tornado.ioloop
import tornado.web
from wtforms_tornado import Form
from wtforms import TextField, StringField, SelectField, FloatField
from tornado.web import RequestHandler
from bokeh.util.browser import view
from bokeh.plotting import figure
from bokeh.layouts import row, column, gridplot, widgetbox, layout
from bokeh.embed import components

path = 'results'
plotkinds = ['tracks','bar','text']
cut_methods = ['default','rank','score']
views = ['binders','promiscuous','summary']

def help_msg():
    msg = '<a>results path not found, enter a folder with your results</a><br>'
    msg += '<a href="%s"> see help page</a>' %wikipage
    return msg

class ControlsForm(Form):
    name = SelectField('name', choices=[])
    path = TextField('path', default='results')
    cutoff = FloatField('cutoff', default=5)
    n = TextField('n', default='2')
    cm = [(i,i) for i in cut_methods]
    cutoff_method = SelectField('cutoff method', choices=cm)
    kinds = [(i,i) for i in plotkinds]
    kind = SelectField('plot kind', choices=kinds)
    views = [(i,i) for i in views]
    view = SelectField('view', choices=views)

class MainHandler(RequestHandler):
    """Handler for main results page"""
    def get(self):
        args = self.request.arguments
        buttons = ''
        self.render('index.html', buttons=buttons)

class GlobalViewHandler(RequestHandler):
    """Handler for showing multiple sequences in a results folder"""

    def get(self):
        args = self.request.arguments
        form = ControlsForm()
        defaultargs = {'path':'results','cutoff':5,'cutoff_method':'rank',
                       'view':'promiscuous','n':2,'kind':'tracks'}
        for k in defaultargs:
            if k in args:
                defaultargs[k] = args[k][0]
        path = defaultargs['path'].strip()
        view = defaultargs['view']
        preds = web.get_predictors(path)

        data = web.get_binder_tables(preds, **defaultargs)
        #add url to prot/seq name
        for k in data:
            data[k] = web.column_to_url(data[k], 'name', '/sequence?path=%s&name=' %path)
        #convert dfs to html
        tables = web.dataframes_to_html(data, classes='tinytable sortable')
        #put tables in tabbed divs
        tables = web.tabbed_html(tables)
        form.path.data = path
        form.cutoff.data = defaultargs['cutoff']
        form.n.data = defaultargs['n']
        form.cutoff_method.data = defaultargs['cutoff_method']
        form.view.data = view

        self.render('global.html', form=form, tables=tables)

class GenomeViewHandler(RequestHandler):
    def get(self):
        args = self.request.arguments
        self.render('genome.html')

class SequenceViewHandler(RequestHandler):
    """Handler for main results page"""
    def get(self):

        args = self.request.arguments
        form = ControlsForm()
        defaultargs = {'path':'results','name':'','cutoff':5,'cutoff_method':'rank',
                       'n':2,'kind':'tracks'}
        for k in defaultargs:
            if k in args:
                defaultargs[k] = args[k][0]
        path = defaultargs['path']
        current_name = defaultargs['name']

        if not os.path.exists(path):
            msg = help_msg()
            self.render('sequence.html', script='', div='', form=form, msg=msg)
            return

        names = web.get_file_lists(path)
        if current_name == '': current_name = names[0]
        form.path.data = path
        form.name.choices = [(i,i) for i in names]
        form.name.data = current_name
        form.cutoff.data = defaultargs['cutoff']
        form.n.data = defaultargs['n']
        form.cutoff_method.data = defaultargs['cutoff_method']
        form.kind.data = defaultargs['kind']

        preds = web.get_predictors(path, current_name)
        plots = web.create_figures(preds, **defaultargs)
        data = web.get_binder_tables(preds, **defaultargs)
        tables = web.dataframes_to_html(data, classes='tinytable sortable')
        tables = web.tabbed_html(tables)
        info = web.dict_to_html(web.get_results_info(preds[0]))

        if len(plots) > 0:
            grid = gridplot(plots, ncols=1, merge_tools=True, sizing_mode='scale_width',
                            toolbar_options=dict(logo=None))
        else:
            grid = ''

        script, div = components(grid)

        self.render('sequence.html', script=script, div=div, form=form, tables=tables,
                    msg='', info=info, name=current_name)

    @staticmethod
    def modify_doc(doc):
        return

settings = dict(
        template_path=os.path.join(os.path.dirname(__file__), "templates"),
        static_path=os.path.join(os.path.dirname(__file__), "static"),
        autoescape=None,
        xsrf_cookies=True,
        debug=True)

def main(port=8888):
    #bokeh_app = Application(FunctionHandler(IndexHandler.modify_doc))
    #bokeh_server = Server({'/main': bokeh_app},
    #                      io_loop=io_loop,
    #                      extra_patterns=[('/', IndexHandler)],
    #                      allow_websocket_origin=['localhost:5006'])
    #bokeh_server.start()
    handlers = [ (r"/", MainHandler),
                 (r"/sequence", SequenceViewHandler),
                 (r"/global", GlobalViewHandler),
                 ]
    app = tornado.web.Application(handlers, **settings)
    #app.listen(8888)
    http_server = tornado.httpserver.HTTPServer(app)
    http_server.listen(port)
    io_loop = tornado.ioloop.IOLoop.current()
    #io_loop.add_callback(view, "http://localhost:8888/")
    view("http://localhost:8888/")
    io_loop.start()


if __name__ == "__main__":
    main()