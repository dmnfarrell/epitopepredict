#!/usr/bin/env python

"""
    epitopepredict server app for viewing results, uses tornado
    Created Sep 2017
    Copyright (C) Damien Farrell
"""

import sys,os,glob
from StringIO import StringIO
import pprint
import pandas as pd
import numpy as np
from epitopepredict import base, web, analysis, config

import tornado.ioloop
import tornado.web
from wtforms_tornado import Form
from wtforms import TextField, StringField, FloatField, IntegerField, BooleanField
from wtforms import SelectField, SelectMultipleField, FileField
from wtforms.validators import DataRequired, Length, ValidationError
from wtforms import widgets
from tornado.web import RequestHandler
from bokeh.util.browser import view
from bokeh.plotting import figure
from bokeh.layouts import row, column, gridplot, widgetbox, layout
from bokeh.embed import components

wikipage = 'https://github.com/dmnfarrell/epitopepredict/wiki/Web-Application'
plotkinds = ['tracks','text','grid']
cut_methods = ['default','rank','score']
views = ['binders','promiscuous','by allele','summary']
opts = config.baseoptions.copy()

def help_msg():
    msg = 'path for results not found, enter an existing folder with your results.  '
    msg += '<a href="%s"> see help page</a>' %wikipage
    return msg

def get_args(args):
    defaults = {'path':'','name':'','cutoff':5,'cutoff_method':'rank', 'pred':'tepitope',
                   'n':2,'kind':'tracks','view':'binders'}
    for k in defaults:
        if k in args:
            defaults[k] = args[k][0]
    return defaults

def str_to_html(s):
    x=''
    for i in s.splitlines():
        x+=i+'<br>'
    return x

def dict_to_html(d):
    x=''
    pp = pprint.PrettyPrinter(depth=2)
    s = pp.pformat(d)
    for i in s.splitlines():
        x+=i+'<br>'
    return x

def is_seqfile(message=u'Wrong format file. Should be fasta or genbank', extensions=None):
    if not extensions:
        extensions = ('fasta', 'faa', 'fa', 'gbk', 'gb')
    def _is_seqfile(form, field):
        if not field.data or field.data.split('.')[-1] not in extensions:
            raise ValidationError(message)
    return _is_seqfile

def exists(message=u'File does not exist'):
    def _exists(form, field):
        if not os.path.exists(field.data):
            raise ValidationError(message)
    return _exists

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
    cached = BooleanField('use cached')

class ConfigForm(Form):
    path = TextField('output path', default='results', validators=[DataRequired()],
                     render_kw={"class": "textbox"})
    pm = [(i,i) for i in base.predictors]
    predictors = SelectMultipleField('predictors', choices=pm,
                                     render_kw={"class": "combobox"})
    mhc1_length = IntegerField('mhc1 length', default=11)
    mhc2_length = IntegerField('mhc2 length', default=15)
    sequence_file = TextField('sequence file',
                              validators=[DataRequired(), is_seqfile(), exists()], default='')
    overwrite = BooleanField('overwrite', default=False, false_values={False, 'n', ''})
    cpus = IntegerField('cpus', default=1)
    ps1 = [(i,i) for i in base.mhc1_presets]
    ps1.insert(0, ('',''))
    mhc1_presets = SelectField('MHC-I presets', choices=ps1, default='')
    ps2 = [(i,i) for i in base.mhc2_presets]
    ps2.insert(0, ('',''))
    mhc2_presets = SelectField('MHC-II presets', choices=ps2, default='')
    p1 = base.get_predictor('iedbmhc1')
    x = [(i,i) for i in p1.getAlleles()]
    mhc1_alleles = SelectMultipleField('MHC-I alleles', choices=x,
                                      render_kw={"class": "combobox"})
    p2 = base.get_predictor('tepitope')
    x = [(i,i) for i in p2.getAlleles()]
    #drballeles = base.getDRBList(mhc2alleles)
    #dqpalleles = base.getDQPList(mhc2alleles)
    mhc2_alleles = SelectMultipleField('MHC-II alleles', choices=x,
                                     render_kw={"class": "combobox"})
    mhc2_alleles.size=5
    iedbmhc1_path = TextField('iedb MHC-I tools path')
    iedbmhc2_path = TextField('iedb MHC-II tools path')

class MainHandler(RequestHandler):
    """Handler for main results page"""
    def get(self):
        args = self.request.arguments
        buttons = ''
        self.render('index.html', buttons=buttons, path='')

class GlobalViewHandler(RequestHandler):
    """Handler for showing multiple sequences in a results folder"""

    def get(self):
        args = self.request.arguments
        form = ControlsForm()
        defaultargs = {'path':'','cutoff':5,'cutoff_method':'rank',
                       'view':'promiscuous','n':2,'cached':1}
        for k in defaultargs:
            if k in args:
                defaultargs[k] = args[k][0]
        path = defaultargs['path'].strip()
        view = defaultargs['view']
        usecached = defaultargs['cached']
        if usecached == 1:
            print ('using cached results')

        if not os.path.exists(path):
            msg = help_msg()
            self.render('global.html', form=form, msg=msg, path=path, status=0)

        preds = web.get_predictors(path)
        data = {}
        if view == 'summary':
            for P in preds:
                if P.data is None: continue
                seqs = web.get_sequences(P)
                #tables = web.sequences_to_html_table(seqs, classes="seqtable")
                pb = P.promiscuousBinders(**defaultargs)
                #print (pb)
                #cl = analysis.find_clusters(b, min_binders=2)
                x = pb.groupby('name').agg({'peptide':np.size,
                                            P.scorekey:np.median}).reset_index()
                x = x.rename(columns={'peptide':'binders'})
                x = x.merge(seqs, on='name', how='right')
                x = web.column_to_url(x, 'name', '/sequence?path=%s&name=' %path)
                data[P.name] = x
        else:
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

        self.render('global.html', form=form, tables=tables, msg='', status=1, path=path)

class GenomeViewHandler(RequestHandler):
    def get(self):
        args = self.request.arguments
        self.render('genome.html')

class SequenceViewHandler(RequestHandler):
    """Handler for main results page"""

    def get(self):
        args = self.request.arguments
        defaultargs = get_args(args)
        form = ControlsForm()
        path = defaultargs['path']
        current_name = defaultargs['name']

        if not os.path.exists(path):
            msg = help_msg()
            self.render('sequence.html', form=form, status=0, name='', msg=msg, path=path)
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
        form.view.data = defaultargs['view']

        preds = web.get_predictors(path, current_name)
        #alleles = web.get_alleles(preds)

        data = web.get_binder_tables(preds, **defaultargs)
        tables = web.dataframes_to_html(data, classes='tinytable sortable')
        tables = web.tabbed_html(tables)
        #info = web.dict_to_html(web.get_results_info(preds[0]))
        info=''
        kind = defaultargs['kind']

        if kind == 'grid':
            seqhtml = web.sequence_to_html_grid(preds, classes="gridtable")
            div = '<div class="scrolled">%s</div>' %seqhtml
            script = ''
        elif kind == 'text':
            seqhtml = web.create_sequence_html(preds, classes="seqtable")
            div = '<div class="scrolled">%s</div>' %seqhtml
            script = ''
        else:
            plots = web.create_figures(preds, **defaultargs)

            if len(plots) > 0:
                grid = gridplot(plots, ncols=1, merge_tools=True, sizing_mode='scale_width',
                                toolbar_options=dict(logo=None))
                script, div = components(grid)
            else:
                script = ''; div = ''

        links = []
        for k in data.keys():
            defaultargs['pred'] = k
            links.append(self.get_url(defaultargs, link=k))

        self.render('sequence.html', script=script, div=div, form=form, tables=tables,
                    msg='', info=info, name=current_name, status=1, links=links, path=path)

    def get_url(self, args, link='download'):
        """Get url from current args"""

        import urllib
        s = '<a href=/download?'
        s += urllib.urlencode(args)
        s += '>%s<a>' %link
        return s

class DownloadHandler(RequestHandler):
    """Download tables as csv"""

    def get(self):
        args = self.request.arguments
        args = get_args(args)
        #args['method'] = 'tepitope'
        filename = args['name']+'_'+args['pred']+'.csv'
        self.set_header ('Content-Type', 'text/csv')
        self.set_header ('Content-Disposition', 'attachment; filename=%s' %filename)
        preds = web.get_predictors(args['path'], args['name'])
        data = web.get_binder_tables(preds, **args)
        out = self.get_csv(data, args['pred'])
        self.write (out)

    def get_csv(self, data, key):
        import io
        if key not in data: return ''
        df=data[key]
        output = io.BytesIO()
        df.to_csv(output, float_format='%.2f')
        csvdata = output.getvalue()
        return csvdata

class MakeConfigHandler(RequestHandler):
    """Make a config file from form"""

    def get(self):
        args = self.request.arguments
        path=''
        for s in opts:
            for i in opts[s]:
                for a in args:
                    if a in opts[s]:
                        opts[s][a] = args[a]

        #print (args)
        o=opts['base']
        if 'overwrite' not in args:
            o['overwrite'] = 'no'
        if 'mhc1_presets' in args and args['mhc1_presets'][0] != '':
            o['mhc1_alleles'] = args['mhc1_presets']
        if 'mhc2_presets' in args and args['mhc2_presets'][0] != '':
            o['mhc2_alleles'] = args['mhc2_presets']
        #get configparser from args to make conf from form
        cp = config.create_config_parser_from_dict(opts)
        out = StringIO()
        cp.write(out)
        conftext = str_to_html(out.getvalue())

        form = ConfigForm(args)
        errors='no errors'
        if form.validate():
            pass
        else:
            errors = 'please fix the following errors:<br>'
            errors += dict_to_html(form.errors)
        self.render('makeconfig.html', form=form, path=path, conftext=conftext, errors=errors)


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
                 (r"/makeconfig", MakeConfigHandler),
                 (r"/download", DownloadHandler)
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