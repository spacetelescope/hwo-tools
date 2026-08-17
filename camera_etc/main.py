import base64
import datetime
import asdf
import yaml
import os
import copy
import numpy as np, astropy.units as u 

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool, Range1d  
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Slider, Div, Select, FileInput
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.layouts import row, column 
from bokeh.io import curdoc

from syotools.spectra.spec_defaults import syn_spectra_library
from syotools.spectra.utils import load_txtfile, load_synfits
from syotools.models import Telescope, Camera, Source, SourcePhotometricExposure

import hdi_help as h 
import synphot as syn
import stsynphot as stsyn

spectra_library = copy.deepcopy(syn_spectra_library)

FLUXUNIT = u.ABmag

hri_source = None
hri_exp = None
hri = None
hwo = None
source1 = ColumnDataSource(data=dict())
source2 = ColumnDataSource(data=dict())
source3 = ColumnDataSource(data=dict())
source4 = ColumnDataSource(data=dict())
suitable_instruments = {}
template_to_start_with = "Flat (AB)"

def update_snr(suitable_instruments, exptime):
    global hri_source
    global hri_exp
    global instrument
    snrs = []
    pivots = []
    names = []
    
    for instrument_name in suitable_instruments:
        bands = suitable_instruments[instrument_name]
        instrument = hwo.instruments[instrument_name]
        bands = sorted(bands, key=lambda x: instrument.configuration["band"][x]["effective_wavelength"])

        snr = []
        pivotwave = []
        instrument.add_exposure(hri_exp)
        hri_exp.source = hri_source
        hri_exp.exptime = exptime
        for band_name in bands:
            try:
                print("Band", band_name)
                hri_exp.calculate_snr(custom_band=band_name)
                band = instrument.configuration["band"][band_name]
                pivotwave.append(band["effective_wavelength"].value)
                snr.append(hri_exp.snr[0].value)
            except (syn.exceptions.DisjointError, syn.exceptions.SynphotError):
                print("Disjoint")
                continue
        snrs.append(np.asarray(snr))
        pivots.append(np.asarray(pivotwave))
        names.append(bands)

    return snrs, pivots, names

def flatten(arr):
    # https://realpython.com/python-flatten-list/
    return [item for row in arr for item in row]

def initialize_setup():
    global hri_source
    global hri_exp
    global hwo

    global source1
    global source2
    global source3
    global source4
    global suitable_instruments


    hwo = Telescope() 
    hwo.set_from_hwome('EAC5')
    suitable_instruments, suitable_bands = hwo.find_instrument_with("filter")

    hri_source = Source() 
    hri_source.set_sed(template_to_start_with, 30., 0., 0.)

    hri_exp = SourcePhotometricExposure() 
    hri_exp.source = hri_source
    hri_exp.verbose = True 

    snrs, pivots, names = update_snr(suitable_instruments, 1 * u.h)

    source1 = ColumnDataSource(data=dict(x=pivots[0], y=snrs[0], desc=names[0]))
    source2 = ColumnDataSource(data=dict(x=pivots[1], y=snrs[1], desc=names[1]))
    source3 = ColumnDataSource(data=dict(x=pivots[2], y=snrs[2], desc=names[2]))
    source4 = ColumnDataSource(data=dict(x=pivots[3], y=snrs[3], desc=names[3]))

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
              x_range=[1200, 23000], y_range=[0, 7], border_fill_color='black', toolbar_location='right')
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
snr_plot.line('x', 'y', source=source4, line_width=3, line_color='green', line_alpha=1.0)
snr_plot.scatter('x', 'y', source=source4, fill_color='white', line_color='green', size=8) 

#hri_source = spectra_library[template_to_start_with]
flux_converted = syn.units.convert_flux(hri_source.sed.waveset, hri_source.sed(hri_source.sed.waveset), FLUXUNIT)

spectrum_template = ColumnDataSource(data=dict(w=hri_source.sed.waveset.value, 
                                               f=flux_converted.value))

sed_plot = figure(height=400, width=700,tools="crosshair,pan,reset,save,box_zoom,wheel_zoom",
              x_range=[800, 24000], y_range=[35, 21], border_fill_color='black', toolbar_location='right')
sed_plot.x_range = Range1d(800, 24000, bounds=(800, 24000)) 
sed_plot.yaxis.axis_label = 'AB Magnitude'
sed_plot.xaxis.axis_label = 'Wavelength [Angstrom]'
sed_plot.line('w','f',line_color='orange', line_width=3, source=spectrum_template, line_alpha=1.0)  



def update_data(attrname, old, new):


    print("You have chosen template ", template.value)

    hwo.effective_diameter = aperture.value * u.m 

    hri_source.set_sed(template.value, magnitude.value, 0., 0., library=spectra_library)

    normalization_band = stsyn.band(spectra_library[template.value].band)
    hri_source.sed.normalize(magnitude.value * u.ABmag, normalization_band) 
    print('Renorming to ', magnitude.value) 
    print('SED Waveunits: ', hri_source.sed.waveset.unit)
    print('SED Fluxunits: ', hri_source.sed(hri_source.sed.waveset).unit)

    flux_converted = syn.units.convert_flux(hri_source.sed.waveset, hri_source.sed(hri_source.sed.waveset), FLUXUNIT)

    spectrum_template.data = {'w':hri_source.sed.waveset.value, 'f':flux_converted.value}    

    snrs, pivots, names = update_snr(suitable_instruments, [exptime.value] * u.hr) 

    source1.data = dict(x=pivots[0], y=snrs[0], desc=names[0]) 
    source2.data = dict(x=pivots[1], y=snrs[1], desc=names[1]) 
    source3.data = dict(x=pivots[2], y=snrs[2], desc=names[2])
    source4.data = dict(x=pivots[3], y=snrs[3], desc=names[3])

    snr_plot.y_range.start = 0
    snr_plot.y_range.end = 1.3*np.max([np.max(flatten(snrs)),5.]) 

    sed_plot.y_range.start = np.max(flux_converted.value)+5.
    sed_plot.y_range.end = np.min(flux_converted.value)-5.

    text = 'Normalized to ' + str(magnitude.value) + ' in the ' + str(spectra_library[template.value].band) + ' band'
    sed_plot.title.text = text
    warning.text = ""

    return source1, source2, source3

# fake source for managing callbacks 
source = ColumnDataSource(data=dict(value=[]))
source.on_change('data', update_data)

aperture = Slider(title="Aperture (meters)", value=hwo.effective_diameter.value, start=4., end=15.0, step=0.1, tags=[4,5,6,6], width=250) 
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


template = Select(title="Template Spectrum", value="Flat (AB)", options=list(spectra_library.keys()), width=250) 

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


for w in [template]:
    w.on_change('value', update_data)


controls = column(children=[aperture, exptime, magnitude, template, upload, warning], sizing_mode='fixed', max_width=300, width=300, height=600) 
controls_tab = TabPanel(child=controls, title='Controls')
info_tab = TabPanel(child=Div(text = h.help()), title='Info')
inputs = Tabs(tabs=[ controls_tab, info_tab], width=300) 
plots = Tabs(tabs=[ TabPanel(child=snr_plot, title='SNR'), TabPanel(child=sed_plot, title='SED') ]) 


curdoc().add_root(row(children=[inputs, plots])) 
curdoc().add_root(source) 
curdoc().add_root(source1) 
