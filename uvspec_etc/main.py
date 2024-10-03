from __future__ import print_function
import numpy as np
import copy 
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool, Range1d 
from bokeh.layouts import row, column
from bokeh.models.widgets import Slider, TextInput, Select, Div 
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.io import curdoc
from bokeh.models.callbacks import CustomJS
import pysynphot as S 
import astropy.constants as const
import get_lumos_spectra
import Telescope as T 
import lumos_help as h 

luvoir = T.Telescope(15., 280., 500.) # set up LUVOIR with 15 meters, T = 280, and diff limit at 500 nm 
lumos = T.Spectrograph() # set up LUVOIR with 10 meters, T = 280, and diff limit at 500 nm 
lumos.set_mode('G120M') 

def simulate_exposure(telescope, spectrograph, wave, flux, exptime): 
    print("Attempting to create an exposure for Telescope: ", telescope.name, telescope.aperture, ' m') 
    print("                                 and Spectrograph: ", spectrograph.name, " in mode ", spectrograph.mode_name) 

    # obtain the interpolated effective areas for the input spectrum 
    aeff_interp = np.interp(wave, spectrograph.wave, spectrograph.aeff, left=0., right=0.) * (telescope.aperture/15.)**2 
    bef_interp = np.interp(wave, spectrograph.wave, spectrograph.bef, left=0., right=0.) # background to use 
    phot_energy = const.h.to('erg s').value * const.c.to('cm/s').value / (wave * 1e-8) # now convert from erg cm^-2 s^-1 A^-1  
    source_counts = flux / phot_energy * aeff_interp * (exptime*3600.) * (wave / lumos.R) 

    source_counts[(wave < spectrograph.lambda_range[0])] = 0. 
    source_counts[(wave > spectrograph.lambda_range[1])] = 0. 

    background_counts = bef_interp / phot_energy * aeff_interp * (exptime*3600.) * (wave / lumos.R) 
    signal_to_noise = source_counts / (source_counts + background_counts)** 0.5 
    return signal_to_noise 

##### START FOR NEW WAY TO GET TEMPLATE SPECTRA 
spec_dict = get_lumos_spectra.add_spectrum_to_library() 
template_to_start_with = 'QSO' 
spec_dict[template_to_start_with].wave 
spec_dict[template_to_start_with].flux # <---- these are the variables you need 

signal_to_noise = simulate_exposure(luvoir, lumos, spec_dict[template_to_start_with].wave, spec_dict[template_to_start_with].flux, 1.0) 

flux_cut = spec_dict[template_to_start_with].flux 
flux_cut[spec_dict[template_to_start_with].wave < lumos.lambda_range[0]] = -999.  
flux_cut[spec_dict[template_to_start_with].wave > lumos.lambda_range[0]] = -999.  

spectrum_template = ColumnDataSource(data=dict(w=spec_dict[template_to_start_with].wave, f=spec_dict[template_to_start_with].flux, \
                                   w0=spec_dict[template_to_start_with].wave, f0=spec_dict[template_to_start_with].flux, \
                                   flux_cut=flux_cut, sn=signal_to_noise)) 

instrument_info = ColumnDataSource(data=dict(wave=lumos.wave, bef=lumos.bef))

# set up the flux plot 
flux_plot = figure(height=400, width=800, 
              tools="crosshair,hover,pan,reset,save,box_zoom,wheel_zoom", outline_line_color='black', 
              x_range=[900, 2000], y_range=[0, 4e-16], toolbar_location='right') 
flux_plot.x_range=Range1d(900,3000,bounds=(900,3000))
flux_plot.y_range=Range1d(0,4e-16,bounds=(0,None)) 
flux_plot.yaxis.axis_label = 'Flux [erg / s / cm^2 / Ang]' 
flux_plot.xaxis.axis_label = 'Wavelength [Angstrom]' 
flux_plot.line('w', 'f', source=spectrum_template, line_width=3, line_color='firebrick', line_alpha=0.7, legend_label='Source Flux')
flux_plot.line('wave', 'bef', source=instrument_info, line_width=3, line_color='darksalmon', line_alpha=0.7, legend_label='Background')

# set up the flux plot 
sn_plot = figure(height=400, width=800, 
              tools="crosshair,hover,pan,reset,save,box_zoom,wheel_zoom", outline_line_color='black', 
              x_range=[900, 2000], y_range=[0, 40], toolbar_location='right')
sn_plot.x_range=Range1d(900,3000,bounds=(900,3000))
sn_plot.y_range=Range1d(0,40,bounds=(0,None)) 
sn_plot.line('w', 'sn', source=spectrum_template, line_width=3, line_color='orange', line_alpha=0.7, legend_label='S/N per resel')
sn_plot.xaxis.axis_label = 'Wavelength [Angstrom]' 
sn_plot.yaxis.axis_label = 'S/N per resel' 

def update_data(attrname, old, new): # use this one for updating pysynphot templates 
   
    print("You have chosen template ", template.value, np.size(spec_dict[template.value].wave)) 
    print('Selected grating = ', grating.value) 
    luvoir.aperture = aperture.value 
    print('Your telescope is set to', luvoir.aperture) 
    lumos.set_mode(grating.value) 

    new_w0 = spec_dict[template.value].wave 
    new_f0 = spec_dict[template.value].flux 

    if ('Blackbody' in template.value):      #<---- update the blackbody curve here. 
       bb = S.BlackBody(bb_temperature.value) 
       bb.convert('flam')
       bb_renorm = bb.renorm(magnitude.value, 'abmag', S.ObsBandpass('galex,fuv')) 
       new_w0 = bb_renorm.wave 
       new_f0 = bb_renorm.flux 
 
    #OOPS, SHOULD USE PYSYNPHOT FOR REDSHIFT HERE, THE NORMALIZATION IS NOT QUITE CORRECT 
    new_w = np.array(new_w0) * (1. + redshift.value)
    new_f = np.array(new_f0) * 10.**( (21.-magnitude.value) / 2.5)
    new_sn = np.nan_to_num(simulate_exposure(luvoir, lumos, new_w, new_f, exptime.value)) 

    flux_cut = copy.deepcopy(new_f) 
    flux_cut[new_w < lumos.lambda_range[0]] = -999.  
    flux_cut[new_w > lumos.lambda_range[1]] = -999.  
    print('RANGE', lumos.lambda_range[0], lumos.lambda_range[1]) 

    new_dict = {'w':new_w, 'f':new_f, 'w0':new_w0, 'f0':new_f0, 'flux_cut':flux_cut, 'sn':new_sn} 
    spectrum_template.data = new_dict 


    # set the axes to autoscale appropriately 
    flux_plot.y_range.start = 0 
    flux_plot.y_range.end = 1.5*np.max(flux_cut)
    sn_plot.y_range.start = 0 
    sn_plot.y_range.end = 1.3*np.max(spectrum_template.data['sn'])
    print('MAX MAX', np.max(spectrum_template.data['f']), np.max(flux_cut)) 

    #instrument_info.data['wave'] = lumos.wave 
    #instrument_info.data['bef'] = lumos.bef  

# fake source for managing callbacks 
source = ColumnDataSource(data=dict(value=[]))
source.on_change('data', update_data)

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
