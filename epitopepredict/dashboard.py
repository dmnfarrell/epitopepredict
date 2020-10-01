#!/usr/bin/env python

"""
    Dashboard app with Bokeh/Panel for epitopepredict
    Created Sep 2019
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

import os
import numpy as np
import pandas as pd
from epitopepredict import plotting, base, peptutils, web, config
from epitopepredict import __version__
#from IPython.display import HTML
import panel as pn
import panel.widgets as pnw
helppage = 'http://epitopepredict.readthedocs.io/en/latest/webapp.html'
css_file=(os.path.join(base.module_path,'static/custom.css'))
js_files = {'sortable': os.path.join(base.module_path,'static/sorttable.js')}
css=''.join(open(css_file,'r').readlines())
pn.extension(raw_css=[css])
pn.config.js_files=js_files

def predictions_dashboard(path):
    """Dashboard for viewing results from epitopepredict runs."""

    #folder_input = pn.widgets.TextInput(name='path', value='../zaire_test', width=400,width_policy='fit')
    #reload_btn = pn.widgets.Button(name='reload',width=100,button_type='primary')

    names = web.get_file_lists(path)
    if names is None:
        return
    preds = web.get_predictors(path,name=names[0])
    print (preds)
    seqname = pnw.Select(name='name', value=names[0], options=names)
    cutoff_slider = pnw.FloatSlider(name='cutoff', value=.95,start=.75,end=.99,step=0.01)
    cutoff_method = pnw.Select(name='cutoff method', value='default', options=['default','rank'])
    n_select = pnw.FloatSlider(name='n',value=1,start=1,end=8,step=1)
    plot_select = pnw.Select(name='plot view', value='tracks', options=['tracks', 'sequence'])
    table_select = pnw.Select(name='table view', value='promiscuous', options=['promiscuous','binders'])
    colorseq_box = pnw.Checkbox(name='color sequences', value=False)

    header = pn.pane.Markdown('__total sequences: %s__' %len(names), css_classes=['main'])
    tables = pn.Tabs(width=900)
    plot = pn.pane.Bokeh(width=800)
    debug = pn.pane.Markdown('test',style={'font-size': '10pt','background-color':'yellow'})
    summary_plot = pn.pane.Bokeh()
    summary_table_tabs = pn.Tabs()
    recalc_button = pnw.Button(name='recalculate',width=200)

    def update_banner():
        """Update the banner"""

        fullpath = os.path.abspath(path)
        banner = pn.Row(pn.pane.Markdown('<h4>epitopepredict: %s</h4> [help](%s) version %s' %(fullpath,helppage,__version__),
                                         css_classes=['divheader'],
                                         sizing_mode='stretch_width'))
        return banner

    def update_header(target, event):
        names = web.get_file_lists(event.new)
        target.object = "_total sequences: %s_" %str(len(names))
        return

    def callback_getpath(event):
        path = os.path.getcwd()
        folder.value = path

    def update_plot(preds,name,cutoff,n,kind):
        """Plot data view"""

        if kind == 'tracks':
            p = plotting.bokeh_plot_tracks(preds,name=name,cutoff=cutoff,n=n,width=1000,title=name)
            plot.object = p
        elif kind == 'sequence':
            p = plotting.bokeh_plot_sequence(preds,name=name,cutoff=cutoff,n=n,width=1000,
                                             title=name,color_sequence=colorseq_box.value)
            plot.object = p
        return p

    def update_tables(preds,name,n):
        """Tabular views of results"""

        P = preds[0]
        view = table_select.value
        tables.clear()
        for P in preds:
            if view == 'promiscuous':
                df = P.promiscuous_binders(n=n,name=name)
            else:
                df = P.get_binders(name=name)
            res = df.to_html(classes="tinytable sortable")
            div = '<div class="scrollingArea">%s</div>' %res
            tables.append((P.name, div))
            #tables.append((P.name,pn.pane.HTML('<p>hddsadsadsasda</p>',width=700)))
        return

    def update(event):
        """Update all elements"""

        name = seqname.value
        n = n_select.value
        cutoff = cutoff_slider.value
        kind = plot_select.value
        debug.object = name
        preds = web.get_predictors(path,name=name)
        update_plot(preds,name=name,cutoff=cutoff,n=n,kind=kind)
        update_tables(preds,name,n)
        return

    def update_summary(path):
        """Summary info for folder"""

        data = web.get_summary_tables(path)
        df = pd.concat(data, sort=True).reset_index()
        #plot = plotting.bokeh_summary_plot(df)
        #summary_plot.object = plot
        summary_table_tabs.clear()
        a = web.aggregate_summary(data)
        div = web.get_scrollable_table(a)
        summary_table_tabs.append(('all', div))
        names = list(data.keys())
        for n in names:
            df = data[n]
            res = df.to_html(classes="tinytable sortable")
            div = '<div class="scrollingArea">%s</div>' %res
            summary_table_tabs.append((n, div))
        return

    @pn.depends(seqname.param.value,n_select.param.value)
    def download_link(name,n):
        if preds is None:
            return
        df = preds[0].promiscuous_binders(n=n,name=name)
        df.to_csv()
        return pn.Pane(HTML('<a>download</a>'),width=700)

    info = pn.pane.Markdown(web.get_readme())
    banner = update_banner()
    update_summary(path)
    #reload_btn.param.watch(load_predictors, 'clicks')
    #reload_btn.param.trigger()
    seqname.param.watch(update, 'value')
    cutoff_slider.param.watch(update, 'value')
    n_select.param.watch(update, 'value')
    table_select.param.watch(update, 'value')
    plot_select.param.watch(update, 'value')
    seqname.param.trigger('options', 'value')

    top = pn.Row(header)#,download_link)
    left = pn.Column(plot,tables,margin=10,sizing_mode='stretch_width')
    right = pn.Column(seqname,cutoff_slider,cutoff_method,n_select,plot_select,
                      table_select,colorseq_box,css_classes=['widget-box'],width=200)
    center = pn.Row(left,right)
    #bottom = pn.Row(table)
    main = pn.Column(top,center)
    summarypane = pn.Column(recalc_button,(pn.Row(summary_table_tabs)))
    tabs = pn.Tabs(('summary', summarypane),('sequence',main),('about',info))
    #tabs.append()
    app = pn.Column(banner,tabs, sizing_mode='stretch_width')
    return app

def run_server(path, port):

    if path == None:
        print ('provide a folder')
        return
    app = predictions_dashboard(path)
    from bokeh.server.server import Server
    def modify_doc(doc):
        return app.server_doc(doc=doc, title='epitopepredict: %s' %path)

    print('Opening application on http://localhost:%s/' %port)
    server = Server({'/': modify_doc}, port=port)
    server.start()
    server.show('/')
    server.run_until_shutdown()
    return

if __name__ == '__main__':
    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-i", "--results", dest="results",
                        help="Results folder", metavar="FILE")
    parser.add_option("-x", "--port", dest="port", default=8000,
                        help="Port for web app, default 8000")

    opts, remainder = parser.parse_args()
    if opts.results == None:
        print ('please provide a results folder')
        sys.exit()
    run_server(opts.results, opts.port)
