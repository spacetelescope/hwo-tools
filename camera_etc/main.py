''' A dead simple ETC for HDST 
'''
from __future__ import print_function
import numpy as np
from bokeh.io import output_file

import phot_compute_snr as phot_etc 

from bokeh.plotting import Figure
from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.models import ColumnDataSource, HoverTool, Paragraph, Range1d 
from bokeh.models.callbacks import CustomJS
from bokeh.layouts import column, row, WidgetBox 
from bokeh.models.widgets import Slider, Tabs, Div, Panel, Select 
from bokeh.io import curdoc
from bokeh.embed import file_html

import Telescope as T 
import hdi_help as h 
import get_hdi_seds 
import pysynphot as S 

luvoir = T.Telescope(10., 280., 500.) # set up LUVOIR with 10 meters, T = 280, and diff limit at 500 nm 
hdi = T.Camera()                     # and HDI camera with default bandpasses 
hdi.set_pixel_sizes(luvoir) 

spec_dict = get_hdi_seds.add_spectrum_to_library() 
template_to_start_with = 'Flat (AB)' 
spec_dict[template_to_start_with].wave 
spec_dict[template_to_start_with].flux # <---- these are the variables you need to make the plot 
spec_dict[template_to_start_with].convert('abmag') 

# set up ColumnDataSources for main SNR plot 
snr = phot_etc.compute_snr(luvoir, hdi, 1., 32.)
source1 = ColumnDataSource(data=dict(x=hdi.pivotwave[2:-3], y=snr[2:-3], desc=hdi.bandnames[2:-3] ))
source2 = ColumnDataSource(data=dict(x=hdi.pivotwave[0:2], y=snr[0:2], desc=hdi.bandnames[0:2]))
source3 = ColumnDataSource(data=dict(x=hdi.pivotwave[-3:], y=snr[-3:], desc=hdi.bandnames[-3:]))

hover = HoverTool(point_policy="snap_to_data", 
        tooltips="""
        <div>
            <div>
                <span style="font-size: 17px; font-weight: bold; color: #696">@desc band</span>
            </div>
            <div>
                <span style="font-size: 15px; font-weight: bold; color: #696">S/N = </span>
                <span style="font-size: 15px; font-weight: bold; color: #696;">@y</span>
            </div>
        </div>
        """
    )


# Set up plot
snr_plot = Figure(plot_height=400, plot_width=800, 
              tools="crosshair,pan,reset,save,box_zoom,wheel_zoom",
              x_range=[120, 2300], y_range=[0, 40], toolbar_location='right')
snr_plot.x_range = Range1d(100, 2300, bounds=(120, 2300)) 
snr_plot.add_tools(hover)
snr_plot.background_fill_color = "beige"
snr_plot.background_fill_alpha = 0.5
snr_plot.yaxis.axis_label = 'SNR'
snr_plot.xaxis.axis_label = 'Wavelength (nm)'
snr_plot.text(5500, 20, text=['V'], text_align='center', text_color='red')

snr_plot.line('x', 'y', source=source1, line_width=3, line_alpha=1.0) 
snr_plot.circle('x', 'y', source=source1, fill_color='white', line_color='blue', size=10)
    
snr_plot.line('x', 'y', source=source2, line_width=3, line_color='orange', line_alpha=1.0)
snr_plot.circle('x', 'y', source=source2, fill_color='white', line_color='orange', size=8) 
    
snr_plot.line('x', 'y', source=source3, line_width=3, line_color='red', line_alpha=1.0)
snr_plot.circle('x', 'y', source=source3, fill_color='white', line_color='red', size=8) 

spectrum_template = ColumnDataSource(data=dict(w=spec_dict[template_to_start_with].wave, f=spec_dict[template_to_start_with].flux, \
                                   w0=spec_dict[template_to_start_with].wave, f0=spec_dict[template_to_start_with].flux))

sed_plot = Figure(plot_height=400, plot_width=800, 
              tools="crosshair,pan,reset,save,box_zoom,wheel_zoom",
              x_range=[120, 2300], y_range=[35, 21], toolbar_location='right')
sed_plot.x_range = Range1d(100, 2300, bounds=(120, 2300)) 
sed_plot.background_fill_color = "beige"
sed_plot.background_fill_alpha = 0.5
sed_plot.yaxis.axis_label = 'AB Mag'
sed_plot.xaxis.axis_label = 'Wavelength (nm)'
sed_plot.line('w','f',line_color='orange', line_width=3, source=spectrum_template, line_alpha=1.0)  




def update_data(attrname, old, new):

    print("You have chosen template ", template.value, np.size(spec_dict[template.value].wave)) 

    luvoir.aperture = aperture.value 
    hdi.set_pixel_sizes(luvoir) # adaptively set the pixel sizes 

    #CHANGE THE SED TEMPLATE IN THE SED PLOT 
    spectrum = spec_dict[template.value]
    band = S.ObsBandpass('johnson,v')
    band.convert('nm') 
    print(spectrum.waveunits.name) 
    ss = spectrum.renorm(magnitude.value+2.5, 'abmag', band) #### OH MY GOD WHAT A HACK!!!! 
    print('Renorming to ', magnitude.value) 
    new_w0 = ss.wave 
    new_f0 = ss.flux 
    new_w = np.array(new_w0) 
    new_f = np.array(new_f0) 
    spectrum_template.data = {'w':new_w, 'f':new_f, 'w0':new_w0, 'f0':new_f0} 

    interp_mags = spec_dict[template.value]
    mag_arr = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]) 
    for pwave,index in zip(hdi.pivotwave,np.arange(10)): 
       mag_arr[index] = ss.sample(pwave) 
       print(index, pwave, ss.sample(pwave), mag_arr[index]) 
    
    snr = phot_etc.compute_snr(luvoir, hdi, exptime.value, mag_arr)
    
    print('SNR') 
    print(mag_arr) 
    print(snr) 
    print('SNR') 

    wave = hdi.pivotwave 
    source1.data = dict(x=wave[2:-3], y=snr[2:-3], desc=hdi.bandnames[2:-3]) 
    source2.data = dict(x=hdi.pivotwave[0:2], y=snr[0:2], desc=hdi.bandnames[0:2]) 
    source3.data = dict(x=hdi.pivotwave[-3:], y=snr[-3:], desc=hdi.bandnames[-3:]) 

    snr_plot.y_range.start = 0
    snr_plot.y_range.end = 1.3*np.max([np.max(snr),5.]) 

    sed_plot.y_range.start = np.min(ss.flux)+5. 
    sed_plot.y_range.end = np.min(ss.flux)-5. 


def fake_function(attrname, old, new):
    print("THIS IS JUST A FAKE CALLBACK FUNCTION", magnitude.value) 

# fake source for managing callbacks 
source = ColumnDataSource(data=dict(value=[]))
source.on_change('data', update_data)

# Set up widgets
aperture= Slider(title="Aperture (meters)", value=6., start=4., end=10.0, step=0.1, tags=[4,5,6,6]) 
aperture_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
aperture.js_on_change('value_throttled', aperture_callback) 

exptime = Slider(title="Exptime (hours)", value=1., start=0.1, end=10.0, step=0.1) 
exptime_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
exptime.js_on_change("value_throttled", exptime_callback) 

magnitude = Slider(title="V Magnitude (AB)", value=28.0, start=20.0, end=35.) 
magnitude_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
magnitude.js_on_change("value_throttled", magnitude_callback) 

template = Select(title="Template Spectrum", value="Flat (AB)", 
                options=["Flat (AB)", "Blackbody (5000K)", "O5V Star", \
                         "B5V Star", "G2V Star", "M2V Star", "Orion Nebula", "Elliptical Galaxy", "Sbc Galaxy", \
                         "Starburst Galaxy", "NGC 1068"]) 
#redshift = Slider(title="Redshift", value=0., start=0.0, end=1.0, step=0.1, callback_policy='mouseup')
#redshift.callback = CustomJS(args=dict(source=source), code="""
#    source.data = { value: [cb_obj.value] }
#""")

for w in [template]: # iterate on changes to parameters 
    w.on_change('value', update_data)

controls = WidgetBox(children=[aperture, exptime, magnitude, template ]) 
controls_tab = Panel(child=controls, title='Controls')
help_tab = Panel(child=Div(text = h.help()), title='Info')
inputs = Tabs(tabs=[ controls_tab, help_tab]) 

plots = Tabs(tabs=[ Panel(child=snr_plot, title='SNR'), Panel(child=sed_plot, title='SED')]) 

# Set up layouts and add to document
curdoc().add_root(row(children=[inputs, plots])) 
curdoc().add_root(source) 
curdoc().add_root(source1) 
