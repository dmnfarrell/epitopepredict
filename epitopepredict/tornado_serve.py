#!/usr/bin/env python

"""
    epitopepredict server app for viewing results, uses tornado
    Created Sep 2017
    Copyright (C) Damien Farrell
"""

import sys,os,glob
from collections import OrderedDict
try:
    from StringIO import StringIO
except:
    from io import StringIO
import pprint
import pandas as pd
import numpy as np
from epitopepredict import base, web, analysis, plotting, config

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

helppage = 'http://epitopepredict.readthedocs.io/en/latest/webapp.html'
plotkinds = ['tracks','text','grid']
cut_methods = ['default','rank','score']
opts = config.baseoptions.copy()

def help_msg():
    msg = 'path for results not found, enter an existing folder with your results.  '
    msg += '<a href="%s"> see help page</a>' %helppage
    return msg

def get_args(args, defaults={'savepath':'','name':'','cutoff':.95,'cutoff_method':'default', 'pred':'tepitope',
                   'n':2,'kind':'tracks','view':'binders'}):
    for k in defaults:
        if k in args:
            defaults[k] = args[k][0].decode("utf-8")
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

class SummaryForm(Form):
    #views = ['summary','promiscuous']
    #name = SelectField('name', choices=[])
    savepath = TextField('savepath', default='results')
    #cutoff = FloatField('cutoff', default=.95)
    #n = TextField('n', default='2')
    #cm = [(i,i) for i in cut_methods]
    #cutoff_method = SelectField('cutoff method', choices=cm)
    #kinds = [(i,i) for i in plotkinds]
    #kind = SelectField('plot kind', choices=kinds)
    #views = [(i,i) for i in views]
    #view = SelectField('table view', choices=views)
    deletecached = BooleanField('delete cached')

class ControlsForm(Form):
    views = ['binders','promiscuous']
    name = SelectField('name', choices=[])
    savepath = TextField('savepath', default='results')
    cutoff = FloatField('cutoff', default=.95)
    n = TextField('n', default='2')
    cm = [(i,i) for i in cut_methods]
    cutoff_method = SelectField('cutoff method', choices=cm)
    kinds = [(i,i) for i in plotkinds]
    kind = SelectField('plot kind', choices=kinds)
    views = [(i,i) for i in views]
    view = SelectField('table view', choices=views)
    #deletecached = BooleanField('delete cached')

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
    x = [(i,i) for i in p1.get_alleles()]
    mhc1_alleles = SelectMultipleField('MHC-I alleles', choices=x,
                                      render_kw={"class": "combobox"})
    p2 = base.get_predictor('tepitope')
    x = [(i,i) for i in p2.get_alleles()]
    #drballeles = base.getDRBList(mhc2alleles)
    #dqpalleles = base.getDQPList(mhc2alleles)
    mhc2_alleles = SelectMultipleField('MHC-II alleles', choices=x,
                                     render_kw={"class": "combobox"})
    mhc2_alleles.size=5
    iedbmhc1_path = TextField('iedb MHC-I tools path')
    iedbmhc2_path = TextField('iedb MHC-II tools path')

class NeoForm(Form):
    savepath = TextField('savepath', default='', validators=[DataRequired()],
                     render_kw={"class": "textbox"})
    #views = ['all','promiscuous','by protein']
    sample = SelectField('sample', choices=[])
    views = ['final','combined']
    views = [(i,i) for i in views]
    view = SelectField('view', choices=views)

class MainHandler(RequestHandler):
    """Handler for main results page"""
    def get(self):
        args = self.request.arguments
        buttons = ''
        self.render('index.html', buttons=buttons, savepath='')

class GlobalViewHandler(RequestHandler):
    """Handler for showing multiple sequences in a results folder"""

    def get(self):
        args = self.request.arguments
        form = SummaryForm()
        defaultargs = {'savepath':'','deletecached':1}
        for k in defaultargs:
            if k in args:
                defaultargs[k] = args[k][0].decode("utf-8")
        savepath = defaultargs['savepath'].strip()
        deletecached = defaultargs['deletecached']
        #if usecached == 1:
        #    print ('using cached results')

        if not os.path.exists(savepath):
            msg = help_msg()
            self.render('global.html', form=form, msg=msg, savepath=savepath, status=0)

        data = web.get_summary_tables(savepath, **defaultargs)

        df = pd.concat(data).reset_index()
        plot = plotting.bokeh_summary_plot(df)
        plots = [plot]

        if len(plots) > 0:
            grid = gridplot(plots, ncols=1, merge_tools=True, sizing_mode='scale_width',
                            toolbar_options=dict(logo=None))
            script, div = components(grid)
        else:
            script = ''; div = ''

        for k in data:
            data[k] = web.column_to_url(data[k], 'name', '/sequence?savepath=%s&name=' %savepath)
        #convert dfs to html
        tables = web.dataframes_to_html(data, classes='tinytable sortable')
        #put tables in tabbed divs
        tables = web.tabbed_html(tables)

        form.savepath.data = savepath
        self.render('global.html', form=form, tables=tables, msg='', script=script, div=div,
                    status=1, savepath=savepath)

class SequenceViewHandler(RequestHandler):
    """Handler for main results page"""

    def get(self):
        args = self.request.arguments
        defaultargs = get_args(args)
        form = ControlsForm()
        savepath = defaultargs['savepath']
        current_name = defaultargs['name']

        if not os.path.exists(savepath):
            msg = help_msg()
            self.render('sequence.html', form=form, status=0, name='', msg=msg, savepath=savepath)
            return

        names = web.get_file_lists(savepath)
        if current_name == '': current_name = names[0]
        form.savepath.data = savepath
        form.name.choices = [(i,i) for i in names]
        form.name.data = current_name
        form.cutoff.data = defaultargs['cutoff']
        form.n.data = defaultargs['n']
        form.cutoff_method.data = defaultargs['cutoff_method']
        form.kind.data = defaultargs['kind']
        form.view.data = defaultargs['view']

        data = web.get_results_tables(path=savepath, **defaultargs)
        tables = web.dataframes_to_html(data, classes='tinytable sortable')
        tables = web.tabbed_html(tables)
        #info = web.dict_to_html(web.get_results_info(preds[0]))
        info=''
        kind = defaultargs['kind']

        preds = web.get_predictors(path=savepath, name=current_name)
        print (preds)

        #if kind == 'grid':
        #    seqhtml = web.sequence_to_html_grid(preds, classes="gridtable")
        #    div = '<div class="scrolled">%s</div>' %seqhtml
        #    script = ''
        if kind == 'text':
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
                    msg='', info=info, name=current_name, status=1, links=links, savepath=savepath)

    def get_url(self, args, link='download'):
        """Get url from current args"""

        import urllib
        s = '<a href=/download?'
        try:
            s += urllib.parse.urlencode(args)
        except:
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

class NeoEpitopeHandler(RequestHandler):
    "Handler for showing neoepitope results"

    def get(self):
        args = self.request.arguments
        form = NeoForm()
        args = get_args(args, defaults={'savepath':'','view':'final','sample':None})
        savepath = args['savepath']
        view = args['view']
        sample = args['sample']
        if not os.path.exists(savepath):
            self.no_data_render(form, msg='no such path')
            return

        labels = pd.read_csv(os.path.join(savepath, 'sample_labels.csv'),index_col=0)
        samples = list(labels.index)
        if sample == None:
            sample = samples[0]

        print (view)
        #get results data
        data = OrderedDict()
        plots = []

        preds = []
        s=sample
        if sample == 'all':
            s=samples[0]
        for p in base.predictors:
            f = os.path.join(savepath, 'results_%s_%s.csv' %(s,p))
            if os.path.exists(f):
                preds.append(p)

        if sample == 'all':
            #get all results together
            for p in preds:
                data['summary'] = labels
                res=[]
                for sample in samples:
                    fname = os.path.join(savepath, 'binders_%s_%s.csv' %(sample,p))
                    b = pd.read_csv(fname,index_col=0)
                    res.append(b)

                res = pd.concat(res).reset_index()
                x = pd.pivot_table(res, index=['name','peptide'], columns=['label'], values='score')
                data[p] = x

        else:
            #get table for variants
            variant_file = os.path.join(savepath, 'variant_effects_%s.csv' %sample)
            if not os.path.exists(variant_file):
                self.no_data_render(form, msg='no such file %s' %variant_file)
                return
            variants = pd.read_csv(variant_file, index_col=0)
            #data['variants'] = variants

            for p in preds:
                binder_file = fname = os.path.join(savepath, 'binders_%s_%s.csv' %(sample,p))
                if not os.path.exists(binder_file):
                    resfile = os.path.join(savepath, 'results_%s_%s.csv' %(sample,p))
                    res = pd.read_csv(resfile,index_col=0)
                    #get binders here if new cutoff used
                else:
                    pb = pd.read_csv(binder_file)
                    pb = pb.drop('label',1)

                #top genes with most peptides?
                t = pb['name'].value_counts()[:20]
                bar = plotting.bokeh_vbar(t, title='top genes with peptide binders: %s' %p, color='#41b6c4')
                plots.append(bar)
                self.add_links(pb)

                #limit table size!
                data[p] = pb[:1000]

            x = variants.variant_class.value_counts()
            pie = plotting.bokeh_pie_chart(x, radius=.25, title='variant classes', height=150)
            plots.append(pie)
            #test = plotting.bokeh_test(height=200)
            v=variants
            x=v.chr.value_counts().sort_index()
            bar2 = plotting.bokeh_vbar(x, title='variants per chromosome', color='#225ea8')
            plots.append(bar2)

        tables = web.dataframes_to_html(data, classes='sortable')
        tables = web.tabbed_html(tables)

        if len(plots) > 0:
            grid = gridplot(plots, ncols=1, nrows=3, merge_tools=True, sizing_mode='scale_width',
                            toolbar_options=dict(logo=None))
            script, div = components(grid)
        else:
            script = ''; div = ''

        form.savepath.data = savepath
        samples.append('all')
        form.sample.choices = [(i,i) for i in samples]
        form.sample.data = sample
        self.render('neoepitope.html', status=1, savepath=savepath, links='', form=form,
                    script=script, div=div, tables=tables)

    def no_data_render(self, form, msg):
        self.render('neoepitope.html', form=form, status=0, name='', msg=msg, savepath='')
        return

    def add_links(self, df):
        web.column_to_url(df, 'name',
                          'http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=')
        web.column_to_url(df, 'transcript_id',
                          'https://www.ensembl.org/id/')

settings = dict(
        template_path=os.path.join(os.path.dirname(__file__), "templates"),
        static_path=os.path.join(os.path.dirname(__file__), "static"),
        autoescape=None,
        xsrf_cookies=True,
        debug=True)

def main(port=8000):
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
                 (r"/download", DownloadHandler),
                 (r"/neoepitope", NeoEpitopeHandler)
                 ]
    app = tornado.web.Application(handlers, **settings)
    #app.listen(8000)
    http_server = tornado.httpserver.HTTPServer(app)
    http_server.listen(port)
    io_loop = tornado.ioloop.IOLoop.current()
    #io_loop.add_callback(view, "http://localhost:8000/")
    view("http://localhost:%s/" %port)
    io_loop.start()


if __name__ == "__main__":
    main()