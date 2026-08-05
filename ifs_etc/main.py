import base64
import datetime
import asdf
import yaml
import os
import copy
import numpy as np, astropy.units as u 
from functools import partial

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d 
from bokeh.layouts import row, column
from bokeh.models.widgets import Slider, Select, Div, Button, FileInput
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.io import curdoc
from bokeh.models.callbacks import CustomJS
import astropy.units as u
import synphot as syn
import stsynphot as stsyn
import ifs_help as h 

from syotools.spectra.spec_defaults import syn_spectra_library
from syotools.spectra.utils import load_txtfile, load_synfits
from syotools.models import Telescope, IFS, Source, SourceIFSExposure

spectra_library = copy.deepcopy(syn_spectra_library)

FLUXUNIT = u.erg / u.s / u.cm**2 / u.AA

hwo = None
instrument = None
ifs_exp = None
snr_results = ColumnDataSource(data={})
spectrum_template = []
instrument_info = ColumnDataSource(data={})
sources = []
suitable_bands = {}
source_num = 0
run_compute = True

hwo = Telescope() 
hwo.set_from_hwome('EAC5')
suitable_instruments, suitable_bands = hwo.find_instrument_with("ifu")

# Set up layouts and add to document
help_text = Div(text = h.help(), width=200) 
help_panel = TabPanel(child=help_text, title='Info') 
source_inputs = Tabs(tabs=[ help_panel], width=300)

grating = Select(title="Grating / Setting", value=list(suitable_bands.keys())[-1], width=200, \
                 options=list(suitable_bands.keys()))

aperture= Slider(title="Aperture (meters)", value=10., start=4., end=15.0, step=0.1, width=200)

exptime = Slider(title="Exposure Time [hr]", value=1.0, start=0.1, end=10.0, step=0.1, width=200)

flux_plot = figure(height=400, width=800, 
              tools="crosshair,hover,pan,reset,save,box_zoom,wheel_zoom", outline_line_color='black', 
              x_range=[900, 2000], y_range=[0, 4e-16], toolbar_location='right') 
flux_plot.x_range=Range1d(900,3000,bounds=(900,3000))
flux_plot.y_range=Range1d(0,4e-16,bounds=(0,None)) 
flux_plot.yaxis.axis_label = 'Flux [erg / s / cm^2 / Ang]' 
flux_plot.xaxis.axis_label = 'Wavelength [Angstrom]' 

sn_plot = figure(height=400, width=800, 
              tools="crosshair,hover,pan,reset,save,box_zoom,wheel_zoom", outline_line_color='black', 
              x_range=[900, 2000], y_range=[0, 40], toolbar_location='right')
sn_plot.x_range=Range1d(900,3000,bounds=(900,3000))
sn_plot.y_range=Range1d(0,40,bounds=(0,None)) 

def update_snr(band_name, instrument_name, exptime):
    global ifs_source
    global ifs_exp
    global instrument
    instrument = hwo.instruments[instrument_name]

    instrument.add_exposure(ifs_exp)
    ifs_exp.exptime = exptime


    ifs_exp.calculate_snr(band=band_name)
    
    snr = ifs_exp.snr[0]
    wave = ifs_exp.wave

    return snr, wave


def update_data(): # use this one for updating synphot templates 
    global sources
    global hwo
    global suitable_bands
    global instrument_info
    # blank out the old list of sources
    ifs_exp.sources = []
    all_fluxes = []
    max_val = np.asarray([0.0])
    for idx, panel in enumerate(sources):
        template = panel.template
        redshift = panel.redshift
        magnitude = panel.magnitude
        bb_temperature = panel.bb_temperature
        print() 
        print() 
        print("You have chosen template ", template.value) 
        print('Selected grating = ', grating.value) 
        print('Your telescope diameter is set to', aperture.value) 
        print('You asked for redshift', redshift.value) 
        hwo.effective_diameter = aperture.value

        ifs_exp.disable()
        
        ifs_source = Source() 
        ifs_source.set_sed(template.value, magnitude.value, redshift.value, 0., library=spectra_library)

        if ('Blackbody' in template.value):      #<---- update the blackbody curve here. 
            wave = np.linspace(100,30000,300) << u.Angstrom
            bb = syn.spectrum.SourceSpectrum(syn.models.BlackBody1D, bb_temperature.value)
            bb.z = redshift.value
            bb = bb.normalize(magnitude.value * u.ABmag, stsyn.band('galex,fuv')) 
            ifs_source.sed = syn.spectrum.SourceSpectrum(syn.models.Empirical1D, points=wave << u.Angstrom, lookup_table=bb(wave))

        flux_converted = syn.units.convert_flux(ifs_source.sed.waveset, ifs_source.sed(ifs_source.sed.waveset), FLUXUNIT)
        print("Flux_converted", flux_converted.value)

        all_fluxes.extend(flux_converted.value)

        spectrum_template[idx].data = dict(w=ifs_source.sed.waveset.value, f=flux_converted.value)

        ifs_exp.add_source(ifs_source)

    max_val = np.ma.max(all_fluxes)

    ifs_exp.verbose = True 

    snr, wave = update_snr(grating.value, suitable_bands[grating.value], exptime.value * u.hr)

    snr_fixed = np.nan_to_num(snr, nan=0)

    snr_results.data = dict(w=wave.value, sn = snr_fixed) 

    instrument_info = ColumnDataSource(data=dict(w=wave.value, bef=instrument.sky(wave).value))

    # set the axes to autoscale appropriately 
    flux_plot.y_range.start = 0 
    flux_plot.y_range.end = 1.5*max_val
    sn_plot.y_range.start = 0 
    sn_plot.y_range.end = 1.3*np.max(snr_results.data['sn'])

    print() 
    print()
    print(instrument)

    snr_compute.label = "--------"
    return snr_results, spectrum_template

def recalculate_snr(event):
    """
    The trick here is that Bokeh only synchronizes the calls at the end of a function, so
    if I want to change the button to "Please wait..." and THEN have it calculate, I have
    to make this call, and have this call fire off another callback with add_next_tick_callback.

    Parameters
    ----------
    newvalues : _type_
        _description_
    """
    global run_compute
    if run_compute:
        snr_compute.label = "Please Wait..."
        curdoc().add_next_tick_callback(partial(update_data))
    run_compute = False

snr_compute = Button(label="Calculate", button_type="primary")
# Can't set up the callback here because we need to define its callback (recalculate_exptime) first.
snr_compute.on_click(recalculate_snr)

def unlock_calcbutton(attrname, old, new):
    global run_compute
    snr_compute.label = "Compute"
    run_compute = True

class source_widget():

    def __init__(self):
        self.template = Select(title="Template Spectrum", value="QSO", options=list(spectra_library.keys()), width=200)
        self.template.on_change("value", unlock_calcbutton)

        self.redshift = Slider(title="Redshift", value=0.0, start=0., end=3.0, step=0.05, width=200)
        self.redshift.on_change("value", unlock_calcbutton)

        self.magnitude = Slider(title="AB Magnitude", value=21., start=15., end=30.0, step=0.1, width=200)
        self.magnitude.on_change("value", unlock_calcbutton)

        self.bb_temperature = Slider(title="Blackbody Temperature [K]", value=10000., start=3000., end=200000.0, step=1000., width=200)
        self.bb_temperature.on_change("value", unlock_calcbutton)

        self.upload = FileInput(accept=[".txt", ".csv", ".ascii", ".fit", ".fits", ".asdf"], title="Upload a Spectrum (.txt or FITS format, 10 MiB max)", directory=False, multiple=False) # 1. list allowed extensions
        self.upload.on_change("filename", self.process_spectrum)
        self.warning = Div(text='<p></p>')

    def process_spectrum(self, attr, old, new):
        global sources
        spectrumhex = self.upload.value
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

            if True:
            #try:
                spectrum = load_synfits({"file": [f"../uploaded/{filename}"], "descs": "uploaded"})

                spectra_library[input_filename] = spectrum
                for source in sources: # add to all the sources' lists
                    if input_filename not in self.template.options:
                        source.template.options.append(input_filename)
                self.template.value = input_filename
                os.remove(f"../uploaded/{filename}") # don't clutter the upload directory
                unlock_calcbutton(None, None, None)
            #except Exception as exc:
                #print(exc)
                #self.warning.text = str(exc)
        else:
            self.warning.text = "File too large"

def add_source_callback(event):
    global source_num
    global sources
    global flux_plot
    # create the new widget
    new_source_widget = source_widget()
    sources.append(new_source_widget)
    source_num += 1
    new_source = column(children=[new_source_widget.template, new_source_widget.redshift, new_source_widget.magnitude, 
                                  new_source_widget.bb_temperature, new_source_widget.upload, new_source_widget.warning], 
                                  sizing_mode='fixed', max_width=300, width=250, height=300)
    print("Column Source", new_source)
    source_panel = TabPanel(child=new_source, title=f"Src {source_num}")

    # add the plot
    spectrum_template.append(ColumnDataSource(data={}))

    for spectrum in spectrum_template:
        flux_plot.line('w', 'f', source=spectrum, line_width=3, line_color='firebrick', line_alpha=0.7, legend_label='Source Flux')
    source_inputs.tabs.append(source_panel)
    unlock_calcbutton(None, None, None)

def initialize_setup():
    global hwo
    global instrument
    global ifs_exp
    global spectrum_template
    global snr_results
    global instrument_info

    ifs_exp = SourceIFSExposure()

    add_source_callback(None)

    update_data()


initialize_setup()

for spectrum in spectrum_template:
    flux_plot.line('w', 'f', source=spectrum, line_width=3, line_color='firebrick', line_alpha=0.7, legend_label='Source Flux')
flux_plot.line('w', 'bef', source=instrument_info, line_width=3, line_color='darksalmon', line_alpha=0.7, legend_label='Background')


sn_plot.line('w', 'sn', source=snr_results, line_width=3, line_color='orange', line_alpha=0.7, legend_label='S/N per resel')
sn_plot.xaxis.axis_label = 'Wavelength [Angstrom]' 
sn_plot.yaxis.axis_label = 'S/N per resel' 

grating.on_change('value', unlock_calcbutton)

# fake source for managing callbacks 
source = ColumnDataSource(data=dict(value=[]))
source.on_change('data', unlock_calcbutton)


add_source = Button(label="Add Source", button_type="success")
add_source.on_click(add_source_callback)

aperture_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
aperture.on_change('value', unlock_calcbutton)

exptime_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
exptime.on_change('value', unlock_calcbutton)
 

exposure_inputs = column(children=[grating, aperture, exptime], sizing_mode='fixed', max_width=300, width=300, height=600 )
exposure_panel = TabPanel(child=exposure_inputs, title='Exposure')
exposure_inputs = Tabs(tabs=[ exposure_panel ], width=300) 

row1 = row(children=[source_inputs, flux_plot])
row2 = row(children=[exposure_inputs, sn_plot])

curdoc().add_root(column(children=[row1, add_source, snr_compute, row2]))
