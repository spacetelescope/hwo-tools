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

source = None
hwo = Telescope()
hwo.set_from_hwome("EAC5")
suitable_instruments, suitable_bands = hwo.find_instrument_with("disperser")

instrument = None
exposure = None
snr_results = ColumnDataSource(data={})
spectrum_template = ColumnDataSource(data={})



FLUXUNIT = u.erg / u.s / u.cm**2 / u.AA

initial_template = "QSO"
initial_magnitude = 21.0 # AB Mag

def update_snr(band_name, instrument_name, exptime):
    global source
    global exposure
    global instrument
    instrument = hwo.instruments[instrument_name]

    instrument.add_exposure(exptime)
    exposure.exptime = exptime

    exposure.calculate_snr(custom_band=band_name)

    # snr is a list, because you can have it run ALL the bands at once
    snr = exposure.snr[0].value
    wave = exposure.wave

    return snr, wave


def initialize_setup():
    global hwo
    global instrument
    global source
    global exposure

    global spectrum_template
    global snr_results
    global instrument_info
    global suitable_bands

    source = Source()
    source.set_sed(initial_template, initial_magnitude, 0., 0.)

    exposure = SourceSpectrographicExposure()
    exposure.source = source
    exposure.verbose = True

    exptime = 1 * u.hr

    initial_band = list(suitable_bands.keys())[-1]
    snr, wave = update_snr(initial_band, suitable_bands[initial_band], exptime)

    # now set up the data for the plots
    # The source SED
    source_flux = syn.units.convert_flux(source.sed.waveset, source.sed(source.sed.waveset), FLUXUNIT)
    spectrum_template = ColumnDataSource(data=dict(w=source.sed.waveset.value, 
                                                   f=source_flux.value)) 

    snr_data = ColumnDataSource(data=dict(w=wave, f=snr))
    background = syn.units.convert_flux(wave, instrument.sky(wave) + exposure.thermal(wave), FLUXUNIT)
    background_data = ColumnDataSource(data=dict(w=wave, f=background.value))
    
