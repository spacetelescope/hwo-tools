import numpy as np, astropy.units as u 

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool, Range1d  
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Slider, Div, Select 
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.layouts import row, column 
from bokeh.io import curdoc

from syotools.spectra.spec_defaults import syn_spectra_library
from syotools.models import Telescope, Camera, Source, SourcePhotometricExposure

import hdi_help as h 
import synphot as syn
import stsynphot as stsyn

hri_source = None
hri_exp = None
hri = None
hwo = None
source1 = ColumnDataSource(data=dict())
source2 = ColumnDataSource(data=dict())
source3 = ColumnDataSource(data=dict())
pivotwave = None
template_to_start_with = "Flat (AB)"

def initialize_setup():
    global hri_source
    global hri_exp
    global hri
    global hwo

    global source1
    global source2
    global source3
    global pivotwave

    hwo = Telescope() 
    hwo.set_from_json('EAC1')
    hri = Camera()   
    hwo.add_camera(hri)

    hri_source = Source() 
    hri_source.set_sed(template_to_start_with, 30., 0., 0.)

    hri_exp = SourcePhotometricExposure() 
    hri_exp.source = hri_source
    hri_exp.verbose = True 
    hri_exp.unknown = 'snr'
    hri.add_exposure(hri_exp) 
    hri_exp._update_snr()

    snr = hri_exp.snr
    pivotwave = np.array(hri.pivotwave[0]) * 10. 

    source1 = ColumnDataSource(data=dict(x=pivotwave[2:-3], y=snr[2:-3], desc=hri.bandnames[2:-3] ))
    source2 = ColumnDataSource(data=dict(x=pivotwave[0:2], y=snr[0:2], desc=hri.bandnames[0:2]))
    source3 = ColumnDataSource(data=dict(x=pivotwave[-3:], y=snr[-3:], desc=hri.bandnames[-3:]))

initialize_setup()





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

snr_plot = figure(height=400, width=700, tools="crosshair,pan,reset,save,box_zoom,wheel_zoom",
              x_range=[1200, 23000], y_range=[0, 10], border_fill_color='black', toolbar_location='right')
snr_plot.x_range = Range1d(1000, 23000, bounds=(1200, 23000)) 
snr_plot.add_tools(hover)
snr_plot.yaxis.axis_label = 'Signal to Noise Ratio'
snr_plot.xaxis.axis_label = 'Wavelength [Angstrom]'

snr_plot.line('x', 'y', source=source1, line_width=3, line_alpha=1.0) 
snr_plot.scatter('x', 'y', source=source1, fill_color='white', line_color='blue', size=10)
snr_plot.line('x', 'y', source=source2, line_width=3, line_color='orange', line_alpha=1.0)
snr_plot.scatter('x', 'y', source=source2, fill_color='white', line_color='orange', size=8) 
snr_plot.line('x', 'y', source=source3, line_width=3, line_color='red', line_alpha=1.0)
snr_plot.scatter('x', 'y', source=source3, fill_color='white', line_color='red', size=8) 

spectrum_template = ColumnDataSource(data=dict(w=syn_spectra_library[template_to_start_with].waveset.value, 
                                               f=syn_spectra_library[template_to_start_with](syn_spectra_library[template_to_start_with].waveset).value)) 

sed_plot = figure(height=400, width=700,tools="crosshair,pan,reset,save,box_zoom,wheel_zoom",
              x_range=[800, 24000], y_range=[35, 21], border_fill_color='black', toolbar_location='right')
sed_plot.x_range = Range1d(800, 24000, bounds=(800, 24000)) 
sed_plot.yaxis.axis_label = 'AB Magnitude'
sed_plot.xaxis.axis_label = 'Wavelength [Angstrom]'
sed_plot.line('w','f',line_color='orange', line_width=3, source=spectrum_template, line_alpha=1.0)  

def update_data(attrname, old, new):

    print("You have chosen template ", template.value) 
    hwo.effective_aperture = aperture.value * u.m 

    hri_source.set_sed(template.value, magnitude.value, 0., 0.)

    normalization_band = stsyn.band(syn_spectra_library[template.value].band)
    hri_source.sed.normalize(magnitude.value * u.ABmag, normalization_band) 
    print('Renorming to ', magnitude.value) 
    print('SED Waveunits: ', hri_source.sed.waveset.unit)
    print('SED Fluxunits: ', hri_source.sed(hri_source.sed.waveset).unit)
    
    spectrum_template.data = {'w':hri_source.sed.waveset.value, 'f':hri_source.sed(hri_source.sed.waveset).value}    

    hri_exp.exptime = [[exptime.value, exptime.value, exptime.value, 
                        exptime.value, exptime.value, exptime.value, 
                        exptime.value, exptime.value, exptime.value, exptime.value], 'hr']
    hri_exp._update_snr() 

    source1.data = dict(x=pivotwave[2:-3], y=hri_exp.snr[2:-3], desc=hri.bandnames[2:-3]) 
    source2.data = dict(x=pivotwave[0:2], y=hri_exp.snr[0:2], desc=hri.bandnames[0:2]) 
    source3.data = dict(x=pivotwave[-3:], y=hri_exp.snr[-3:], desc=hri.bandnames[-3:]) 

    snr_plot.y_range.start = 0
    snr_plot.y_range.end = 1.3*np.max([np.max(hri_exp.snr.value),5.]) 

    sed_plot.y_range.start = np.min(hri_source.sed(hri_source.sed.waveset).value)+5. 
    sed_plot.y_range.end = np.min(hri_source.sed(hri_source.sed.waveset).value)-5. 
    text = 'Normalized to ' + str(magnitude.value) + ' in the ' + str(syn_spectra_library[template.value].band) + ' band'
    sed_plot.title.text = text

    return source1, source2, source3

# fake source for managing callbacks 
source = ColumnDataSource(data=dict(value=[]))
source.on_change('data', update_data)

aperture = Slider(title="Aperture (meters)", value=6., start=4., end=10.0, step=0.1, tags=[4,5,6,6], width=250) 
aperture_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
aperture.js_on_change('value_throttled', aperture_callback) 

exptime = Slider(title="Exptime (hours)", value=1., start=0.1, end=20.0, step=0.1, width=250) 
exptime_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
exptime.js_on_change("value_throttled", exptime_callback) 

magnitude = Slider(title="Magnitude (AB)", value=30.0, start=20.0, end=35., step=0.1, width=250) 
magnitude_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
magnitude.js_on_change("value_throttled", magnitude_callback) 

template = Select(title="Template Spectrum", value="Flat (AB)", options=list(syn_spectra_library.keys()), width=250) 

for w in [template]: 
    w.on_change('value', update_data)

controls = column(children=[aperture, exptime, magnitude, template ], sizing_mode='fixed', max_width=300, width=300, height=600) 
controls_tab = TabPanel(child=controls, title='Controls')
info_tab = TabPanel(child=Div(text = h.help()), title='Info')
inputs = Tabs(tabs=[ controls_tab, info_tab], width=300) 
plots = Tabs(tabs=[ TabPanel(child=snr_plot, title='SNR'), TabPanel(child=sed_plot, title='SED') ]) 

curdoc().add_root(row(children=[inputs, plots])) 
curdoc().add_root(source) 
curdoc().add_root(source1) 
