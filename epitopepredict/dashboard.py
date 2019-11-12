
"""
    dashboard app with panel for epitopepredict
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
from IPython.display import HTML
import panel as pn
import panel.widgets as pnw
helppage = 'http://epitopepredict.readthedocs.io/en/latest/webapp.html'
css_file=(os.path.join(base.module_path,'static/custom.css'))
css=''.join(open(css_file,'r').readlines())
pn.extension(raw_css=[css])

def predictions_dashboard(path):
    """Dashboard for viewing results from epitopepredict runs."""
    
    #folder_input = pn.widgets.TextInput(name='path', value='../zaire_test', width=400,width_policy='fit')
    #reload_btn = pn.widgets.Button(name='reload',width=100,button_type='primary')    
    
    names = web.get_file_lists(path)
    if names is None:
        return
    preds = web.get_predictors(path,name=names[0])
    print (preds)
    seqname = pn.widgets.Select(name='name', value=names[0], options=names)
    cutoff_slider = pn.widgets.FloatSlider(name='cutoff', value=.95,start=.75,end=.99,step=0.01)
    cutoff_method = pn.widgets.Select(name='cutoff method', value='default', options=['default','rank'])
    n_select = pn.widgets.FloatSlider(name='n',value=1,start=1,end=8,step=1)
    plot_select = pn.widgets.Select(name='plot view', value='tracks', options=['tracks', 'text'])
    table_select = pn.widgets.Select(name='table view', value='promiscuous', options=['promiscuous','binders'])
    header = pn.pane.Markdown('__total sequences: %s__' %len(names), css_classes=['main'])
    tables = pn.Tabs(width=900)
    plot = pn.pane.Bokeh(width=800)
    debug = pn.pane.Markdown('test',style={'font-size': '10pt','background-color':'yellow'})
    summary_pane = pn.pane.Bokeh()
    
    def update_banner():
        """Update the banner"""
        banner = pn.Row(pn.pane.Markdown('<h4>epitopepredict</h4> [help](%s) version %s' %(helppage,__version__),
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
            #p=plotting.plot_tracks([P],name=name,cutoff=cutoff,n=n,width=700,height=250)   
            plot.object = p
        else:        
            p = plotting.bokeh_plot_grid(P,name=name,width=700,height=250)            
            #plot.object = pn.Pane(div,width=700)
        return p

    def update_table(preds,name,n):
        """Tabular views of results"""

        P = preds[0]
        view = table_select.value
        tabdata = []
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
        #print (tabdata)        
        #tables.value = tabdata
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
        update_table(preds,name,n)
        return
    
    def update_summary(path):
        """Summary info for folder"""
        
        data = web.get_summary_tables(path)
        df = pd.concat(data, sort=True).reset_index()
        plot = plotting.bokeh_summary_plot(df)   
        summary_pane.object = plot
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
    seqname.param.trigger('options', 'value')

    top = pn.Row(header)#,download_link)
    left = pn.Column(plot,tables,margin=10)
    right = pn.Column(seqname,cutoff_slider,cutoff_method,n_select,plot_select,
                      table_select,css_classes=['widget-box'],width=200)
    center = pn.Row(left,right)
    #bottom = pn.Row(table)
    main = pn.Column(top,center)
    tabs = pn.Tabs(('summary', summary_pane),('sequence',main),('about',info))
    #tabs.append()    
    app = pn.Column(banner,tabs, sizing_mode='stretch_width')
    return app

path = 'zaire_test'
app = predictions_dashboard(path)
app.servable()