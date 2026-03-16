import base64
import os
import datetime
import copy

import yaml
import numpy as np
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d 
from bokeh.layouts import row, column
from bokeh.models.widgets import Slider, Select, Div, FileInput
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.io import curdoc
from bokeh.models.callbacks import CustomJS
import astropy.units as u
import synphot as syn
import stsynphot as stsyn

import uvi_help as h 
from syotools.spectra.spec_defaults import syn_spectra_library
from syotools.spectra.utils import load_txtfile, load_synfits
from syotools.models import Telescope, Spectrograph, Source, SourceSpectrographicExposure

spectra_library = copy.deepcopy(syn_spectra_library)

hwo = None
uvi = None
uvi_exp = None
snr_results = ColumnDataSource(data={})
spectrum_template = ColumnDataSource(data={})
instrument_info = ColumnDataSource(data={})

def initialize_setup():
    global hwo
    global uvi
    global uvi_exp
    global spectrum_template
    global snr_results
    global instrument_info

    hwo = Telescope() 
    hwo.set_from_sei('EAC1')
    uvi = Spectrograph()
    uvi.set_from_sei("UVI")
    hwo.add_spectrograph(uvi)              

    template_to_start_with = 'QSO' 

    uvi_source = Source() 
    uvi_source.set_sed(template_to_start_with, 21., 0., 0.)

    uvi_exp = SourceSpectrographicExposure() 
    uvi_exp.source = uvi_source
    uvi_exp.verbose = True 
    uvi_exp.unknown = 'snr'
    uvi.add_exposure(uvi_exp) 
    uvi_exp._update_snr(uvi_source) 

    spectrum_template = ColumnDataSource(data=dict(w=uvi_source.sed.waveset.value, f=uvi_source.sed(uvi_source.sed.waveset).value)) 
    print(' flux = ', uvi_source.sed(uvi_source.sed.waveset))

    snr_results = ColumnDataSource(data=dict(w=uvi.wave.value, sn = uvi_exp.snr.value)) 

    instrument_info = ColumnDataSource(data=dict(wave=uvi.wave.value, bef=uvi.bef.value))

initialize_setup()

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
    uvi_source.set_sed(template.value, magnitude.value, redshift.value, 0., library=spectra_library)

    uvi_exp.exptime = [[exptime.value, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0], 'hr'] 
    uvi_exp.source = uvi_source
    uvi_exp.verbose = True 
    uvi_exp.unknown = 'snr'
    uvi.add_exposure(uvi_exp) 
    uvi.mode = str.split(grating.value,' ')[0] #<-- becuase text after the space not in mode keys  


    if ('Blackbody' in template.value):      #<---- update the blackbody curve here. 
       wave = np.linspace(100,30000,300) << u.Angstrom
       bb = syn.spectrum.SourceSpectrum(syn.models.BlackBody1D, bb_temperature.value)
       bb.z = redshift.value
       bb = bb.normalize(magnitude.value * u.ABmag, stsyn.band('galex,fuv')) 
       uvi_source.sed = syn.spectrum.SourceSpectrum(syn.models.Empirical1D, points=wave, lookup_table=bb(wave))

    uvi_exp._update_snr(uvi_source) 

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
    print(uvi)
    return snr_results, spectrum_template

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

upload = FileInput(accept=[".txt", ".csv", ".ascii", ".fit", ".fits", ".asdf"], title="Upload a Spectrum (.txt or FITS format, 10 MiB max)", directory=False, multiple=False) # 1. list allowed extensions
warning = Div(text='<p></p>')

def process_spectrum(attr, old, new):
    global template
    global spectra_library
    spectrumhex = upload.value
    if len(spectrumhex) < 13981013: #10 MiB in base64 5. Set a file size limit
        spectrumdata = base64.b64decode(spectrumhex, validate=True)
        keyword = spectrumdata[0:6].decode()
        input_filename = new
        if len(input_filename) > 44:
            input_filename = new[0:44]
        
        filetype = "unknown"
        if keyword == "SIMPLE": # 2. Validate the file type, don't trust Content-Type header
            filetype = "fits"
        elif keyword[0:5] == "#ASDF":
            filetype = "asdf"
        elif keyword[0:5] == "%YAML":
            filetype = "yaml"
        else:
            filetype = "txt"

        filename = f"file_{datetime.datetime.now().isoformat()}.{filetype}" # 3. Change the filename to something generated by the application. 6. store files... outside of the webroot 
        with open(f"../uploaded/{filename}", "wb") as outfile:
            outfile.write(spectrumdata)

        try:
            spectrum = load_synfits({"file": [f"../uploaded/{filename}"], "descs": "uploaded"})

            spectra_library[input_filename] = spectrum
            if input_filename not in template.options:
                template.options.append(input_filename)
            template.value = input_filename
            os.remove(f"../uploaded/{filename}") # don't clutter the upload directory
            update_data("","","")
        except Exception as exc:
            warning.text = str(exc)
    else:
        warning.text = "File too large"

upload.on_change("filename", process_spectrum)


# iterate on changes to parameters 
for w in [template, grating]:  w.on_change('value', update_data)
 
# Set up layouts and add to document
help_text = Div(text = h.help(), width=200) 
source_inputs = column(children=[template, redshift, magnitude, bb_temperature, upload, warning], sizing_mode='fixed', max_width=300, width=250, height=300)
controls_panel = TabPanel(child=source_inputs, title='Source') 
help_panel = TabPanel(child=help_text, title='Info') 
source_inputs = Tabs(tabs=[ controls_panel, help_panel], width=300) 
#source_inputs = Tabs(tabs=[ controls_panel], width=300) 

exposure_inputs = column(children=[grating, aperture, exptime], sizing_mode='fixed', max_width=300, width=300, height=600 )
exposure_panel = TabPanel(child=exposure_inputs, title='Exposure')
exposure_inputs = Tabs(tabs=[ exposure_panel ], width=300) 

row1 = row(children=[source_inputs, flux_plot])
row2 = row(children=[exposure_inputs, sn_plot])

curdoc().add_root(column(children=[row1, row2]))
