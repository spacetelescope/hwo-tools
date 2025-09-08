import numpy as np
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d 
from bokeh.layouts import row, column
from bokeh.models.widgets import Slider, Select, Div 
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.io import curdoc
from bokeh.models.callbacks import CustomJS
import astropy.units as u
import synphot as syn
import stsynphot as stsyn
import uvi_help as h 

from syotools.spectra.spec_defaults import syn_spectra_library
from syotools.models import Telescope, Spectrograph, Source, SourceSpectrographicExposure

hwo = Telescope() 
hwo.set_from_json('EAC1')
uvi = Spectrograph()   
hwo.add_spectrograph(uvi)              

template_to_start_with = 'QSO' 

uvi_source = Source() 
uvi_source.set_sed(template_to_start_with, 21., 0., 0.)

uvi_exp = SourceSpectrographicExposure() 
uvi_exp.source = uvi_source
uvi_exp.verbose = True 
uvi_exp.unknown = 'snr'
uvi.add_exposure(uvi_exp) 
uvi_exp._update_snr() 

spectrum_template = ColumnDataSource(data=dict(w=uvi_source.sed.waveset.value, f=uvi_source.sed(uvi_source.sed.waveset).value)) 
print(' flux = ', uvi_source.sed(uvi_source.sed.waveset))

snr_results = ColumnDataSource(data=dict(w=uvi.wave.value, sn = uvi_exp.snr.value)) 

instrument_info = ColumnDataSource(data=dict(wave=uvi.wave.value, bef=uvi.bef.value))

flux_plot = figure(height=400, width=800, 
              tools="crosshair,hover,pan,reset,save,box_zoom,wheel_zoom", outline_line_color='black', 
              x_range=[900, 2000], y_range=[0, 4e-16], toolbar_location='right') 
flux_plot.x_range=Range1d(900,3000,bounds=(900,3000))
flux_plot.y_range=Range1d(0,4e-16,bounds=(0,None)) 
flux_plot.yaxis.axis_label = 'Flux [erg / s / cm^2 / Ang]' 
flux_plot.xaxis.axis_label = 'Wavelength [Angstrom]' 
flux_plot.line('w', 'f', source=spectrum_template, line_width=3, line_color='firebrick', line_alpha=0.7, legend_label='Source Flux')
flux_plot.line('wave', 'bef', source=instrument_info, line_width=3, line_color='darksalmon', line_alpha=0.7, legend_label='Background')

sn_plot = figure(height=400, width=800, 
              tools="crosshair,hover,pan,reset,save,box_zoom,wheel_zoom", outline_line_color='black', 
              x_range=[900, 2000], y_range=[0, 40], toolbar_location='right')
sn_plot.x_range=Range1d(900,3000,bounds=(900,3000))
sn_plot.y_range=Range1d(0,40,bounds=(0,None)) 
sn_plot.line('w', 'sn', source=snr_results, line_width=3, line_color='orange', line_alpha=0.7, legend_label='S/N per resel')
sn_plot.xaxis.axis_label = 'Wavelength [Angstrom]' 
sn_plot.yaxis.axis_label = 'S/N per resel' 

def update_data(attrname, old, new): # use this one for updating pysynphot templates 
    
    print() 
    print() 
    print("You have chosen template ", template.value) 
    print('Selected grating = ', grating.value) 
    print('Your telescope is set to', aperture.value) 
    print('You asked for redshift', redshift.value) 
    hwo.effective_aperture = aperture.value     
    
    uvi_source = Source() 
    uvi_source.set_sed(template.value, magnitude.value, redshift.value, 0.)

    uvi_exp.exptime = [[exptime.value, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0], 'hr'] 
    uvi_exp.source = uvi_source
    uvi_exp.verbose = True 
    uvi_exp.unknown = 'snr'
    uvi.add_exposure(uvi_exp) 
    uvi.mode = str.split(grating.value,' ')[0] #<-- becuase text after the space not in mode keys  


    if ('Blackbody' in template.value):      #<---- update the blackbody curve here. 
       uvi_source = syn.spectrum.SourceSpectrum(syn.blackbody.BlackBody1D, bb_temperature) 
       uvi_source = uvi_source.normalize(magnitude.value * u.ABmag, stsyn.band('galex,fuv')) 

    uvi_exp._update_snr() 

    snr_fixed = np.nan_to_num(uvi_exp.snr.value, nan=0)

    spectrum_template.data = dict(w=uvi_source.sed.waveset.value, f=uvi_source.sed(uvi_source.sed.waveset).value) 
    snr_results.data = dict(w=uvi.wave.value, sn = snr_fixed) 

    # set the axes to autoscale appropriately 
    flux_plot.y_range.start = 0 
    flux_plot.y_range.end = 1.5*np.max(uvi_source.sed(uvi_source.sed.waveset).value)
    sn_plot.y_range.start = 0 
    sn_plot.y_range.end = 1.3*np.max(snr_results.data['sn'])

    print() 
    print() 

# fake source for managing callbacks 
source = ColumnDataSource(data=dict(value=[]))
source.on_change('data', update_data)

# Set up widgets and their callbacks (faking the mouseup policy via "source" b/c functional callback doesn't do that. 
template = Select(title="Template Spectrum", value="QSO", options=list(syn_spectra_library.keys()), width=200)

redshift = Slider(title="Redshift", value=0.0, start=0., end=3.0, step=0.05, width=200)
redshift_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
redshift.js_on_change("value_throttled", redshift_callback) 

magnitude = Slider(title="AB Magnitude", value=21., start=15., end=30.0, step=0.1, width=200)
magnitude_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
magnitude.js_on_change("value_throttled", magnitude_callback) 

bb_temperature = Slider(title="Blackbody Temperature [K]", value=20000., start=8000., end=200000.0, step=1000., width=200)
temperature_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
bb_temperature.js_on_change("value_throttled", temperature_callback) 

grating = Select(title="Grating / Setting", value="G150M (R = 30,000)", width=200, \
                 options=["G120M (R = 30,000)", "G150M (R = 30,000)", "G180M (R = 30,000)", "G155L (R = 5,000)", "G145LL (R = 500)"])

aperture= Slider(title="Aperture (meters)", value=6., start=4., end=10.0, step=0.1, width=200)
aperture_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
aperture.js_on_change("value_throttled", aperture_callback) 

exptime = Slider(title="Exposure Time [hr]", value=1.0, start=0.1, end=10.0, step=0.1, width=200)
exptime_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
exptime.js_on_change("value_throttled", exptime_callback) 

# iterate on changes to parameters 
for w in [template, grating]:  w.on_change('value', update_data)
 
# Set up layouts and add to document
help_text = Div(text = h.help(), width=200) 
source_inputs = column(children=[template, redshift, magnitude, bb_temperature], sizing_mode='fixed', max_width=300, width=250, height=300)
controls_panel = TabPanel(child=source_inputs, title='Source') 
help_panel = TabPanel(child=help_text, title='Info') 
source_inputs = Tabs(tabs=[ controls_panel, help_panel], width=300) 
source_inputs = Tabs(tabs=[ controls_panel], width=300) 

exposure_inputs = column(children=[grating, aperture, exptime], sizing_mode='fixed', max_width=300, width=300, height=600 )
exposure_panel = TabPanel(child=exposure_inputs, title='Exposure')
exposure_inputs = Tabs(tabs=[ exposure_panel ], width=300) 

row1 = row(children=[source_inputs, flux_plot])
row2 = row(children=[exposure_inputs, sn_plot])

curdoc().add_root(column(children=[row1, row2]))
