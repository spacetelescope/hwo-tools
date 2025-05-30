import os

import numpy as np
from bokeh.plotting import figure
from bokeh.embed import components
from bokeh.models import ColumnDataSource, HoverTool, Paragraph, Range1d  
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Slider, Div, Select 
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.layouts import row, column 
from bokeh.io import curdoc, output_file

import astropy.units as u
import astropy.constants as c
import synphot as syn
from syotools.sci_eng_interface import read_json 

import pyEDITH as pE
import pickle
import eacy
from synphot import SourceSpectrum, SpectralElement, Observation
from synphot.models import BlackBodyNorm1D, Empirical1D
from synphot import units as synphot_u

import catalog

eacy.load_instrument("CI", plotting=True, verbose=True) # load the coronagraph instrument given in the yaml files. This also calculates the throughput, and we can tell it to plot. 
eacy.load_detector("IFS", plotting=True, verbose=True) # load the detector instrument given in the yaml files. This also calculates the throughput, and we can tell it to plot.

param_snr=10

target_catalog = catalog.load_catalog()

parameters = {}

def update_scene(updates):

    parameters.update(updates)
    scene = pE.AstrophysicalScene()
    scene.load_configuration(parameters)
    scene.calculate_zodi_exozodi(parameters)
    scene.validate_configuration()

    return scene

def update_observation(updates):

    parameters.update(updates)

    observation = pE.Observation() # define the observation object
    observation.load_configuration(parameters) # load the specified configuration in the parameters dict 
    observation.set_output_arrays()
    observation.validate_configuration()
    
    return observation

def update_observatory(updates, observation, scene):
    parameters.update(updates)

    observatory_config = pE.parse_input.get_observatory_config(parameters)
    observatory = pE.ObservatoryBuilder.create_observatory(observatory_config)
    pE.ObservatoryBuilder.configure_observatory(
        observatory, parameters, observation, scene
    )
    observatory.validate_configuration()

def print_observatory(observatory):
    print("Parameters in observatory:")
    for key in observatory.__dict__.keys():
        if key in ["telescope", "coronagraph", "detector"]:
            print("-->", key)
            subkeys = observatory.__dict__[key].__dict__.keys()
            for subkey in subkeys:
                print("    -->", subkey)
        else:
            print("-->",key)

def recalculate_snr(observation, scene, observatory):
    pE.calculate_exposure_time_or_snr(observation, scene, observatory, verbose=False)


def load_initial():
    """
    Load initial parameters for the ETC on load
    """
    global parameters

    # observation parameters
    # set up wavelengths
    parameters["wavelength"] = np.linspace(0.3, 1.7, 1000)
    parameters["nlambd"] = len(parameters["wavelength"]) # number of wavelengths
    parameters["snr"] = param_snr*np.ones_like(parameters["wavelength"]) # the SNR you want for each spectral bin 
    parameters["CRb_multiplier"] = 2. # factor to multiply the background by (used for differential imaging)
    #parameters["photometric_aperture_radius"] = None#0.85 # radius of the photometric aperture in units of lambda/D
    parameters["psf_trunc_ratio"] = 0.3 # truncate the off-axis PSFs at this level 

    parameters["regrid_wavelength"] = True # set the flag to do this. We also need to specify a few other parameters.
    parameters["spectral_resolution"] = np.array([50, 140, 50]) #np.array([140])  # we're going to define three spectral channels. These are the spectral resolutions for each channel. i.e. all spectral bins in a given channel will have a fixed resolution.
    parameters["channel_bounds"] = np.array([0.4, 1.0]) #np.array([]) # specify the boundaries between the channels in um

    observation = update_observation({})

    # The Astrophysical Scene
    # STAR
    parameters["Lstar"] = 1. # luminosity of the star in solar luminosities
    parameters["distance"] = 10. # distance to the system in pc
    # Note: we can work in either mag or flux units. Let's choose to work in flux units. 

    parameters["stellar_angular_diameter"] = 1.5e-4 # angular diameter of the star = 0.01 lam/D == about 1.5e-4 arcsec for lam=0.5 um and D~6.5 m
    parameters["Fstar_10pc"] = Fstar_obs_10pc.value # Fstar modeled as a blackbody (see above)
    #parameters["FstarV_10pc"] = v_band_flux.value # V-band flux of the star; you don't have to pre-calculate this, ETC will do it for you

    # PLANET
    parameters["separation"] = 0.1 # planetary separation in arcsec
    parameters["Fp/Fs"] = FpFs # 1e-8 for testing (bright planet)
    parameters["Fp_min/Fs"] = 1e-11 # minimum detectable planet 


    # SCENE
    parameters["nzodis"] = 3. # number of zodis for exozodi estimate
    parameters["ra"] = 176.6292 # approximate ra of HD 102365. WARNING: do not use this number for science. 
    parameters["dec"] = -40.5003 # approximate dec of HD 102365. WARNING: do not use this number for science. 

    scene = update_scene({})

    # Observatory parameters
    parameters["observing_mode"] = "IFS" # ETC should use IFS mode
    parameters["observatory_preset"] = "EAC1" # tells ETC to use EAC1 yaml files throughputs
    parameters["IFS_eff"]  = 1. # extra throughput of the IFS 
    parameters["npix_multiplier"] = np.ones_like(parameters["wavelength"]) # number of detector pixels per spectral bin
    parameters["noisefloor_PPF"] = 30 # post processing factor of 30 is a good realistic value for this

    observatory = update_observatory({}, observation, scene)



