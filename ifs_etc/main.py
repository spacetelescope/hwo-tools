import numpy as np
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d 
from bokeh.layouts import row, column
from bokeh.models.widgets import Slider, Select, Div , Button
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.io import curdoc
from bokeh.models.callbacks import CustomJS
import astropy.units as u
import synphot as syn
import stsynphot as stsyn
import ifs_help as h 

from syotools.spectra.spec_defaults import syn_spectra_library
from syotools.models import Telescope, Spectrograph, Source, SourceIFSExposure

hwo = None
ifs = None
ifs_exp = None
snr_results = ColumnDataSource(data={})
spectrum_template = []
instrument_info = ColumnDataSource(data={})
sources = []
source_num = 0

# Set up layouts and add to document
help_text = Div(text = h.help(), width=200) 
help_panel = TabPanel(child=help_text, title='Info') 
source_inputs = Tabs(tabs=[ help_panel], width=300)

grating = Select(title="Grating / Setting", value="G150M (R = 30,000)", width=200, \
                 options=["G120M (R = 30,000)", "G150M (R = 30,000)", "G180M (R = 30,000)", "G155L (R = 5,000)", "G145LL (R = 500)"])

aperture= Slider(title="Aperture (meters)", value=6., start=4., end=10.0, step=0.1, width=200)

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

def update_data(attrname, old, new): # use this one for updating synphot templates 
    global sources
    # blank out the old list of sources
    ifs_exp.sources = []
    max_val = 0
    for idx, panel in enumerate(sources):
        template = panel["template"]
        redshift = panel["redshift"]
        magnitude = panel["magnitude"]
        bb_temperature = panel["bb_temperature"]
        print() 
        print() 
        print("You have chosen template ", template.value) 
        print('Selected grating = ', grating.value) 
        print('Your telescope is set to', aperture.value) 
        print('You asked for redshift', redshift.value) 
        hwo.effective_aperture = aperture.value  
        
        ifs_source = Source() 
        ifs_source.set_sed(template.value, magnitude.value, redshift.value, 0.)

        max_val = np.max([np.max(ifs_source.sed(ifs_source.sed.waveset).value),max_val])

        if ('Blackbody' in template.value):      #<---- update the blackbody curve here. 
            wave = np.linspace(100,30000,300) << u.Angstrom
            bb = syn.spectrum.SourceSpectrum(syn.models.BlackBody1D, bb_temperature.value)
            bb.z = redshift.value
            bb = bb.normalize(magnitude.value * u.ABmag, stsyn.band('galex,fuv')) 
            ifs_source.sed = syn.spectrum.SourceSpectrum(syn.models.Empirical1D, points=wave, lookup_table=bb(wave))

        spectrum_template[idx].data = dict(w=ifs_source.sed.waveset.value, f=ifs_source.sed(ifs_source.sed.waveset).value)
        print(len(spectrum_template))

        ifs_exp.add_source(ifs_source)

    ifs_exp.exptime = [[exptime.value, exptime.value, exptime.value, exptime.value, exptime.value, exptime.value, exptime.value, exptime.value, exptime.value, exptime.value], 'hr'] 
    ifs_exp.verbose = True 
    ifs_exp.unknown = 'snr'
    ifs.add_exposure(ifs_exp) 
    ifs.mode = str.split(grating.value,' ')[0] #<-- because text after the space not in mode keys  



    ifs_exp.calculate()

    snr_fixed = np.nan_to_num(ifs_exp.snr.value, nan=0)

    snr_results.data = dict(w=ifs.wave.value, sn = snr_fixed) 

    # set the axes to autoscale appropriately 
    flux_plot.y_range.start = 0 
    flux_plot.y_range.end = 1.5*max_val
    sn_plot.y_range.start = 0 
    sn_plot.y_range.end = 1.3*np.max(snr_results.data['sn'])

    print() 
    print()
    print(ifs)
    return snr_results, spectrum_template

def create_source_widget():

    source_widget = {}
    # Set up widgets and their callbacks (faking the mouseup policy via "source" b/c functional callback doesn't do that. 
    source_widget["template"] = Select(title="Template Spectrum", value="QSO", options=list(syn_spectra_library.keys()), width=200)
    source_widget["template"].on_change("value", update_data)

    source_widget["redshift"] = Slider(title="Redshift", value=0.0, start=0., end=3.0, step=0.05, width=200)
    redshift_callback = CustomJS(args=dict(source=source_widget["redshift"]), code="""
        source.data = { value: [cb_obj.value] }
        """)
    source_widget["redshift"].js_on_change("value_throttled", redshift_callback) 

    source_widget["magnitude"] = Slider(title="AB Magnitude", value=21., start=15., end=30.0, step=0.1, width=200)
    magnitude_callback = CustomJS(args=dict(source=source_widget["magnitude"]), code="""
        source.data = { value: [cb_obj.value] }
        """)
    source_widget["magnitude"].js_on_change("value_throttled", magnitude_callback) 

    source_widget["bb_temperature"] = Slider(title="Blackbody Temperature [K]", value=20000., start=8000., end=200000.0, step=1000., width=200)
    temperature_callback = CustomJS(args=dict(source=source_widget["bb_temperature"]), code="""
        source.data = { value: [cb_obj.value] }
        """)
    source_widget["bb_temperature"].js_on_change("value_throttled", temperature_callback)


    return source_widget


def add_source_callback(event):
    global source_num
    global sources
    global flux_plot
    # create the new widget
    source_widget = create_source_widget()
    sources.append(source_widget)
    source_num += 1
    new_source = column(children=[source_widget["template"], source_widget["redshift"], source_widget["magnitude"], source_widget["bb_temperature"]], sizing_mode='fixed', max_width=300, width=250, height=300)
    print("Column Source", new_source)
    source_panel = TabPanel(child=new_source, title=f"Src {source_num}")

    # add the plot
    spectrum_template.append(ColumnDataSource(data={}))

    for spectrum in spectrum_template:
        flux_plot.line('w', 'f', source=spectrum, line_width=3, line_color='firebrick', line_alpha=0.7, legend_label='Source Flux')
    source_inputs.tabs.append(source_panel)
    update_data(None, None, None)

def initialize_setup():
    global hwo
    global ifs
    global ifs_exp
    global spectrum_template
    global snr_results
    global instrument_info

    hwo = Telescope() 
    hwo.set_from_json('EAC1')
    ifs = Spectrograph()   
    hwo.add_spectrograph(ifs)       

    ifs_exp = SourceIFSExposure()        

    add_source_callback(None)

    snr_results = ColumnDataSource(data=dict(w=ifs.wave.value, sn = ifs_exp.snr.value)) 

    instrument_info = ColumnDataSource(data=dict(wave=ifs.wave.value, bef=ifs.bef.value))

initialize_setup()

for spectrum in spectrum_template:
    flux_plot.line('w', 'f', source=spectrum, line_width=3, line_color='firebrick', line_alpha=0.7, legend_label='Source Flux')
flux_plot.line('wave', 'bef', source=instrument_info, line_width=3, line_color='darksalmon', line_alpha=0.7, legend_label='Background')


sn_plot.line('w', 'sn', source=snr_results, line_width=3, line_color='orange', line_alpha=0.7, legend_label='S/N per resel')
sn_plot.xaxis.axis_label = 'Wavelength [Angstrom]' 
sn_plot.yaxis.axis_label = 'S/N per resel' 

grating.on_change('value', update_data)

# fake source for managing callbacks 
source = ColumnDataSource(data=dict(value=[]))
source.on_change('data', update_data)


add_source = Button(label="Add Source", button_type="success")
add_source.on_click(add_source_callback)


aperture_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
aperture.js_on_change("value_throttled", aperture_callback) 

exptime_callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
exptime.js_on_change("value_throttled", exptime_callback)
 

exposure_inputs = column(children=[grating, aperture, exptime], sizing_mode='fixed', max_width=300, width=300, height=600 )
exposure_panel = TabPanel(child=exposure_inputs, title='Exposure')
exposure_inputs = Tabs(tabs=[ exposure_panel ], width=300) 

row1 = row(children=[source_inputs, flux_plot])
row2 = row(children=[exposure_inputs, sn_plot])

curdoc().add_root(column(children=[row1, add_source, row2]))
