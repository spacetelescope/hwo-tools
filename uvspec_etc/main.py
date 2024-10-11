import numpy as np
import copy 
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool, Range1d, Div 
from bokeh.layouts import row, column
from bokeh.models.widgets import Slider, Select
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.io import curdoc
from bokeh.models.callbacks import CustomJS
import pysynphot as S 
import astropy.constants as const
import uvi_help as h 
import astropy.units as u
from syotools.models import Telescope, Spectrograph

# create Telescope, Camera, and Exposure objects 
hwo, uvi = Telescope(), Spectrograph() 
hwo.set_from_json('EAC3')
hwo.add_spectrograph(uvi)
uvi_exp = uvi.create_exposure()
uvi.mode='G150M' 

template_to_start_with = 'QSO' 
uvi_exp.sed_id = template_to_start_with 
uvi_exp.renorm_sed(21. * u.ABmag, bandpass='galex,fuv', waveunits='angstrom', fluxunits='flam')
uvi_exp.exptime[1]['value'] = [1.0*3600., 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]

uvi_exp.enable()
uvi_sed, uvi_snr = uvi_exp.recover('sed', 'snr')

spectrum_template = ColumnDataSource(data=dict(w=uvi_exp.exp_wave, f=uvi_exp.exp_flux, sn=uvi_exp.exp_snr)) 
instrument_info = ColumnDataSource(data=dict(wave=uvi.wave, bef=uvi.bef))
sed_cds = ColumnDataSource(data=dict(w=uvi_exp._sed.wave, f=uvi_exp._sed.flux)) 

# set up the flux plot 
flux_plot = figure(height=400, width=800, 
              tools="crosshair,hover,pan,reset,save,box_zoom,wheel_zoom", outline_line_color='black', 
              x_range=[900, 2000], y_range=[0, 1e-14], toolbar_location='right') 
flux_plot.x_range, flux_plot.y_range = Range1d(900,3000,bounds=(900,3000)), Range1d(0,1e-14,bounds=(0,None)) 
flux_plot.yaxis.axis_label, flux_plot.xaxis.axis_label = 'Flux [erg / s / cm^2 / Ang]', 'Wavelength [Angstrom]' 
flux_plot.line('w', 'f', source=sed_cds, line_width=2, line_color='cyan', line_alpha=0.7, legend_label='Source SED')
flux_plot.line('w', 'f', source=spectrum_template, line_width=3, line_color='firebrick', line_alpha=0.7, legend_label='Observed Flux')
flux_plot.line('wave', 'bef', source=instrument_info, line_width=3, line_color='darksalmon', line_alpha=0.7, legend_label='Background')

# set up the snr plot 
sn_plot = figure(height=400, width=800, 
              tools="crosshair,hover,pan,reset,save,box_zoom,wheel_zoom", outline_line_color='black', 
              x_range=[900, 2000], y_range=[0, 25], toolbar_location='right')
sn_plot.x_range=Range1d(900,3000,bounds=(900,3000))
sn_plot.y_range=Range1d(0,25,bounds=(0,None)) 
sn_plot.line('w', 'sn', source=spectrum_template, line_width=3, line_color='orange', line_alpha=0.7, legend_label='S/N per resel')
sn_plot.xaxis.axis_label = 'Wavelength [Angstrom]' 
sn_plot.yaxis.axis_label = 'S/N per resel' 
sn_plot.yaxis.major_label_overrides = {0:' ', 5:'                5'} 

def update_cds(attrname, old, new): # use this one for updating pysynphot templates 

    print("You have chosen template ", template.value) 
    hwo.aperture = aperture.value
    print('Your telescope is set to', hwo.aperture)
    print('Selected grating = ', grating.value,  uvi.list_of_gratings[grating.value]) 
    uvi.mode = uvi.list_of_gratings[grating.value] 
    print(' ') 

    uvi_sed, uvi_snr = uvi_exp.recover('sed', 'snr')

    uvi_exp.sed_id = template.value 
    uvi_exp.renorm_sed(magnitude.value * u.ABmag, bandpass='galex,fuv', waveunits='angstrom', fluxunits='flam')
    uvi_exp.exptime[1]['value'] = [exptime.value * 3600., 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    
    uvi_sed, uvi_snr = uvi_exp.recover('sed', 'snr')

    #need to enable extinction here also, probablu inside SYOTools 
    
    if ('Blackbody' in template.value):      #<---- update the blackbody curve here. 
       print(" Blackbody not implemented yet") 
       #bb = S.BlackBody(bb_temperature.value) 
       #bb.convert('flam')
       #bb_renorm = bb.renorm(magnitude.value, 'abmag', S.ObsBandpass('galex,fuv')) 
       #spec_with_extinction = bb_renorm * S.Extinction(reddening.value, 'mwavg')  
  
    spectrum_template.data = dict(w=uvi_exp.exp_wave, f=uvi_exp.exp_flux, sn=uvi_exp.exp_snr)
    instrument_info = ColumnDataSource(data=dict(wave=uvi.wave, bef=uvi.bef))
    sed_cds.data = dict(w=uvi_exp._sed.wave, f=uvi_exp._sed.flux) 

    flux_plot.y_range.start = 0   #   set the axes to autoscale appropriately 
    flux_plot.y_range.end = 1.5*np.max(np.array(uvi_exp.exp_flux)) 


# fake source for managing callbacks 
source = ColumnDataSource(data=dict(value=[]))
source.on_change('data', update_cds)

# Set up widgets and their callbacks (faking the mouseup policy via "source" b/c functional callback doesn't do that. 
template = Select(title="Template Spectrum", value="QSO", width=200, \
                options=["Flat in F_lambda", "QSO", "10 Myr Starburst", "O5V Star",\
                         "G2V Star", "G191B2B (WD)", "GD71 (WD)", "GD153 (WD)", \
                         "Classical T Tauri", "M1 Dwarf", "Orion Nebula", \
                         "Starburst, No Dust", "Starburst, E(B-V) = 0.6", \
                         "Galaxy with f_esc, HI=1, HeI=1", "Galaxy with f_esc, HI=0.001, HeI=1",\
			 "Blackbody"]) 

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

reddening = Slider(title="Extinction E(B-V)", value=0., start=0., end=2.0, step=0.05, width=200)
reddening_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
reddening.js_on_change("value_throttled", reddening_callback) 

bb_temperature = Slider(title="Blackbody Temperature [K]", value=20000., start=8000., end=200000.0, step=1000., width=200)
temperature_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
bb_temperature.js_on_change("value_throttled", temperature_callback) 

grating = Select(title="Grating / Setting", value="G150M (R = 30,000)", width=200, \
                 options=["G120M (R = 30,000)", "G150M (R = 30,000)", "G180M (R = 30,000)", "G155L (R = 5,000)", "G145LL (R = 500)"])

aperture= Slider(title="Aperture (meters)", value=8., start=4., end=10.0, step=0.1, width=200)
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
for w in [template, grating]:  w.on_change('value', update_cds)
 
# Set up the layouts and add to document
source_controls = column(children=[template, redshift, magnitude, reddening, bb_temperature], sizing_mode='fixed', max_width=300, width=300, height=300)
controls_panel = TabPanel(child=source_controls, title='Source') 
help_text = Div(text = h.help(), width=100, height=300, max_width=300, sizing_mode='inherit') 
help_panel = TabPanel(child=help_text, title='Info') 
#source_inputs = Tabs(tabs=[ controls_panel, help_panel], width=300) 
source_inputs = Tabs(tabs=[ controls_panel], width=300) 

exposure_inputs = column(children=[grating, aperture, exptime], sizing_mode='fixed', max_width=300, width=300, height=600 )
exposure_panel = TabPanel(child=exposure_inputs, title='Exposure')
exposure_inputs = Tabs(tabs=[ exposure_panel], width=300) 

row1 = row(children=[source_inputs, flux_plot])
row2 = row(children=[exposure_inputs, sn_plot])

curdoc().add_root(column(children=[row1, row2]))
