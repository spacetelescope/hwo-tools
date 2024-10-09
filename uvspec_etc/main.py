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

# create a Telescope, Camera, and Exposure 
hwo, uvi = Telescope(), Spectrograph() 
hwo.set_from_json('EAC3')
hwo.add_spectrograph(uvi)
uvi_exp = uvi.create_exposure()
uvi.mode='G120M' 

#spec_dict = get_uvi_spectra.add_spectrum_to_library() # obtain a dictionary with the pysynphot spectra in it  
template_to_start_with = 'QSO' 
uvi_exp.sed_id = template_to_start_with 
uvi_exp.renorm_sed(21. * u.ABmag, bandpass='galex,fuv', waveunits='angstrom', fluxunits='flam')
uvi_exp.exptime[1]['value'] = [1.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]

uvi_exp.verbose = True
hwo.verbose = True
snr = uvi_exp.recover('snr')
uvi_exp.enable()
uvi_sed, uvi_snr = uvi_exp.recover('sed', 'snr')

spectrum_template = ColumnDataSource(data=dict(w=uvi_exp.exp_wave, f=uvi_exp.exp_flux, sn=uvi_exp.exp_snr)) 
instrument_info = ColumnDataSource(data=dict(wave=uvi.wave, bef=uvi.bef))

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

# set up the snr plot 
sn_plot = figure(height=400, width=800, 
              tools="crosshair,hover,pan,reset,save,box_zoom,wheel_zoom", outline_line_color='black', 
              x_range=[900, 2000], y_range=[0, 40], toolbar_location='right')
sn_plot.x_range=Range1d(900,3000,bounds=(900,3000))
sn_plot.y_range=Range1d(0,40,bounds=(0,None)) 
sn_plot.line('w', 'sn', source=spectrum_template, line_width=3, line_color='orange', line_alpha=0.7, legend_label='S/N per resel')
sn_plot.xaxis.axis_label = 'Wavelength [Angstrom]' 
sn_plot.yaxis.axis_label = 'S/N per resel' 
sn_plot.yaxis.major_label_overrides = {0:' ', 5:'                5'} 


def update_cds(attrname, old, new): # use this one for updating pysynphot templates 

    print("You have chosen template ", template.value) 
    print('Selected grating = ', grating.value)
    hwo.aperture = aperture.value
    print('Your telescope is set to', hwo.aperture)
    uvi.mode = grating.value

    uvi_exp.sed_id = template.value 
    uvi_exp.renorm_sed(magnitude.value * u.ABmag, bandpass='galex,fuv', waveunits='angstrom', fluxunits='flam')
    uvi_exp.exptime[1]['value'] = [exptime.value, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    
    snr = uvi_exp.recover('snr')
    uvi_exp.enable()
    uvi_sed, uvi_snr = uvi_exp.recover('sed', 'snr')
    
    spectrum_template.data = dict(w=uvi_exp.exp_wave, f=uvi_exp.exp_flux, sn=uvi_exp.exp_snr)
    instrument_info = ColumnDataSource(data=dict(wave=uvi.wave, bef=uvi.bef))

    #set the axes to autoscale appropriately 
    flux_plot.y_range.start = 0 
    flux_plot.y_range.end = 1.5*np.max(np.array(uvi_exp.exp_flux)) 
    sn_plot.y_range.start = 0
    #sn_plot.y_range.end = 1.3*np.max(np.array(uvi_exp.exp_snr)) 

def update_data(attrname, old, new): # use this one for updating pysynphot templates 
   
    print("You have chosen template ", template.value, np.size(spec_dict[template.value].wave)) 
    print('Selected grating = ', grating.value) 
    hwo.aperture = aperture.value 
    print('Your telescope is set to', hwo.aperture) 
    uvi.mode = grating.value  

    print('Applying reddening for E(B-V)') 
    spec_with_extinction = spec_dict[template.value] * S.Extinction(reddening.value, 'mwavg')  

    if ('Blackbody' in template.value):      #<---- update the blackbody curve here. 
       bb = S.BlackBody(bb_temperature.value) 
       bb.convert('flam')
       bb_renorm = bb.renorm(21., 'abmag', S.ObsBandpass('galex,fuv')) 
       spec_with_extinction = bb_renorm * S.Extinction(reddening.value, 'mwavg')  
 
    spec_with_extinction_renorm = spec_with_extinction.renorm(21., 'abmag', S.ObsBandpass('galex,fuv')) 

    new_w0 = spec_with_extinction_renorm.wave 
    new_f0 = spec_with_extinction_renorm.flux 

    #OOPS, SHOULD USE PYSYNPHOT FOR REDSHIFT HERE, THE NORMALIZATION IS NOT QUITE CORRECT 
    new_w = np.array(new_w0) * (1. + redshift.value)
    new_f = np.array(new_f0) * 10.**( (21.-magnitude.value) / 2.5)
    new_sn = np.nan_to_num(simulate_exposure(hwo, uvi, new_w, new_f, exptime.value)) 

    flux_cut = copy.deepcopy(new_f) 
    flux_cut[new_w < uvi.wrange[0]] = -999.  
    flux_cut[new_w > uvi.wrange[1]] = -999.  
    print('RANGE', uvi.wrange[0], uvi.wrange[1]) 

    new_dict = {'w':new_w, 'f':new_f, 'w0':new_w0, 'f0':new_f0, 'flux_cut':flux_cut, 'sn':new_sn} 
    spectrum_template.data = new_dict 

    # set the axes to autoscale appropriately 
    flux_plot.y_range.start = 0 
    flux_plot.y_range.end = 1.5*np.max(flux_cut)
    sn_plot.y_range.start = 0 
    sn_plot.y_range.end = 1.3*np.max(spectrum_template.data['sn'])
    print('MAX MAX', np.max(spectrum_template.data['f']), np.max(flux_cut)) 

# fake source for managing callbacks 
source = ColumnDataSource(data=dict(value=[]))
#source.on_change('data', update_data)
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
#for w in [template, grating]:  w.on_change('value', update_data)
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
