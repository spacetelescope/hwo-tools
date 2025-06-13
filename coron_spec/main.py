import os

import numpy as np
import scipy as sc
from bokeh.plotting import figure
from bokeh.embed import components
from bokeh.models import ColumnDataSource, Paragraph, Range1d, RadioButtonGroup 
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Slider, Div, Select 
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.layouts import row, column, layout
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

#eacy.load_instrument("CI", plotting=True, verbose=True) # load the coronagraph instrument given in the yaml files. This also calculates the throughput, and we can tell it to plot. 
#eacy.load_detector("IFS", plotting=True, verbose=True) # load the detector instrument given in the yaml files. This also calculates the throughput, and we can tell it to plot.

''' 
Notes to fix:

1. Because each changed value needs to be escaped from the list it's in, we can and should easily divide up which parts of PyEDITH need to be rerun
2. Implement a "Calculate" button because PyEdith takes time, instead of any change to data.
3. Implement the load_star and load_planet features
4. Get the hover-text working
'''

param_snr=10
EACS = ["EAC1", "EAC2", "EAC3"]

target_planet, target_star = catalog.load_catalog()

parameters = {}
scene = None
observation = None

def compute_blackbody_photon_flux(temp, wavelengths, dist):
    """Generate photon flux density (photon/s/cm^2/um) for a blackbody at 1 cm^2."""
    bb = SourceSpectrum(BlackBodyNorm1D, temperature=temp)
    

    flux_photlam = syn.units.convert_flux(wavelengths, bb(wavelengths), u.photon / (u.s * u.cm**2 * u.nm)) #.value * 1/(u.s * u.cm**2 * u.AA) # photons/s/cm2/A
    #flux_photlam = flux_photlam.to(1 / (u.s * u.cm**2 * u.nm)) # convert to 1/nm
    flux_photlam = flux_photlam * (1000*u.pc/dist)**2

    return flux_photlam

def update_scene():
    global parameters
    #parameters.update(updates)

    scene = pE.AstrophysicalScene()
    scene.load_configuration(parameters)
    scene.calculate_zodi_exozodi(parameters)
    scene.validate_configuration()

    print("Planetflux", scene.Fp_over_Fs)

    return scene

def update_observation():
    global parameters
    #parameters.update(updates)

    observation = pE.Observation() # define the observation object
    observation.load_configuration(parameters) # load the specified configuration in the parameters dict 
    observation.set_output_arrays()
    observation.validate_configuration()
    
    return observation

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

# Ordinarily, these would be separate, but at the moment all changes here would seem to affect observatory, observation, and scene
# it is, particularly, unclear what 
def update_calculation(updates):
    global parameters

    if "startemplate" in updates:
        parameters = load_star(updates)
    if "planettemplate" in updates:
        parameters = load_planet(updates)
    if "newsnr" in updates:
        parameters["snr"] = updates["newsnr"][0]

    if "eacnum" in parameters:
        parameters["observatory_preset"] = EACS[updates["eacnum"][0]]
    else:
        parameters["observatory_preset"] = "EAC1" # tells ETC to use EAC1 yaml files throughputs

    scene = update_scene()

    observation = update_observation()

    if parameters["regrid_wavelength"] is True:
        scene.regrid_spectra(parameters, observation)


    observatory_config = pE.parse_input.get_observatory_config(parameters)

    observatory = pE.ObservatoryBuilder.create_observatory(observatory_config)
    pE.ObservatoryBuilder.configure_observatory(
        observatory, parameters, observation, scene
    )
    observatory.validate_configuration()
    print_observatory(observatory)

    return observatory, scene, observation

def recalculate_exptime(updates, old, new):
    print(updates)
    print(old)
    print(new)
    observatory, scene, observation = update_calculation(new)

    pE.calculate_exposure_time_or_snr(observation, scene, observatory, verbose=True)

    good = np.where(observation.exptime < 1e7 * u.s)

    return observation.wavelength[good], observation.exptime[good]

def recalculate_snr(updates):
    observatory, scene, observation = update_calculation(updates)

    pE.calculate_exposure_time_or_snr(observation, scene, observatory, mode="signal_to_noise", verbose=False)

    good = np.where(~np.isinf(observation.fullsnr))

    return observation.wavelength, observation.fullsnr

def load_initial():
    """
    Load initial parameters for the ETC on load
    """
    global parameters

    # observation parameters
    # set up wavelengths
    parameters["wavelength"] = np.linspace(0.3, 1.7, 1000)
    parameters["nlambd"] = len(parameters["wavelength"]) # number of wavelengths
    parameters["snr"] = param_snr * np.ones_like(parameters["wavelength"]) # the SNR you want for each spectral bin 
    parameters["CRb_multiplier"] = 2. # factor to multiply the background by (used for differential imaging)
    #parameters["photometric_aperture_radius"] = None#0.85 # radius of the photometric aperture in units of lambda/D
    parameters["psf_trunc_ratio"] = 0.3 # truncate the off-axis PSFs at this level 

    parameters["regrid_wavelength"] = True # set the flag to do this. We also need to specify a few other parameters.
    parameters["spectral_resolution"] = np.array([70, 140, 70]) #np.array([140])  # we're going to define three spectral channels. These are the spectral resolutions for each channel. i.e. all spectral bins in a given channel will have a fixed resolution.
    parameters["channel_bounds"] = np.array([0.4, 1.0]) #np.array([]) # specify the boundaries between the channels in um


    # The Astrophysical 
    # STAR
    parameters["Lstar"] = 1. # luminosity of the star in solar luminosities
    parameters["distance"] = 10. # distance to the system in pc
    # Note: we can work in either mag or flux units. Let's choose to work in flux units. 
    Fstar_obs_10pc = compute_blackbody_photon_flux(5770, parameters["wavelength"] << u.um, parameters["distance"] << u.pc)

    parameters["stellar_angular_diameter"] = 1.5e-4 # angular diameter of the star = 0.01 lam/D == about 1.5e-4 arcsec for lam=0.5 um and D~6.5 m
    parameters["Fstar_10pc"] = Fstar_obs_10pc.value # Fstar modeled as a blackbody (see above)
    #parameters["FstarV_10pc"] = v_band_flux.value # V-band flux of the star; you don't have to pre-calculate this, ETC will do it for you
    print("Star", parameters["Fstar_10pc"])

    # PLANET
    # Most of these input spectra are albedo for a system. Actual conversion would require knowing the system separation and relative sizes of the star and planet. OR we make that an input parameter.
    parameters["delta_mag"] = 15 * u.ABmag
    FpFs = target_planet["Earth"]["spectrum"](parameters["wavelength"] << u.um)
    # normalize
    contrast = 1e-8
    FpFs = FpFs/np.max(FpFs)
    Fplan_obs = Fstar_obs_10pc*contrast
    FpFs = (Fplan_obs / Fstar_obs_10pc).value

    from matplotlib import pyplot as plt
    plt.plot(parameters["wavelength"] * u.um, Fstar_obs_10pc)
    plt.plot(parameters["wavelength"] * u.um, Fplan_obs)
    plt.yscale("log")

    plt.show()

    #FpFs = FpFs * np.median(Fstar_obs_10pc) - syn.units.convert_flux(parameters["wavelength"], FpFs, parameters["delta_mag"]
    print("Planet", FpFs, Fplan_obs)
    #flux_planet = Fstar_obs_10pc * FpFs * parameters["delta_mag"]
    #parameters["F0"] = flux_planet.value
    parameters["separation"] = 0.1 # planetary separation in arcsec
    parameters["Fp/Fs"] = FpFs # 1e-8 for testing (bright planet)
    parameters["Fp_min/Fs"] = 1e-11 # minimum detectable planet 
    print("Star and Planet", parameters["Fstar_10pc"], parameters["Fp/Fs"])

    # SCENE
    parameters["nzodis"] = 3. # number of zodis for exozodi estimate
    parameters["ra"] = 176.6292 # approximate ra of HD 102365. WARNING: do not use this number for science. 
    parameters["dec"] = -40.5003 # approximate dec of HD 102365. WARNING: do not use this number for science. 

    # Observatory parameters
    parameters["observing_mode"] = "IFS" # ETC should use IFS mode
    if "eacnum" in parameters:
        parameters["observatory_preset"] = EACS[parameters["eacnum"]]
    else:
        parameters["observatory_preset"] = "EAC1" # tells ETC to use EAC1 yaml files throughputs
    parameters["IFS_eff"]  = 1. # extra throughput of the IFS 
    parameters["npix_multiplier"] = np.ones_like(parameters["wavelength"]) # number of detector pixels per spectral bin
    parameters["noisefloor_PPF"] = 30 # post processing factor of 30 is a good realistic value for this

    wavelength, exptime = recalculate_exptime(parameters, [], [])

    return wavelength, exptime

def load_planet(sourceID, star):
    global parameters
    reflect_planet = target_planet[sourceID]["spectrum"]

    flux_planet = star * reflect_planet(parameters["wavelengths"])
    parameters["F0"] = flux_planet

    return flux_planet

def load_star(sourceID):
    global parameters
    flux_star = target_star[sourceID]["spectrum"]

    parameters["Fstar_10pc"] = flux_star(parameters["wavelengths"])

    print("Star Flux", flux_star)

    return flux_star


wavelength, exptime = load_initial()
obsdata = ColumnDataSource(data={"wavelength": wavelength, "exptime": exptime})

print("Resulting Exposure Time", exptime)

inputs = ColumnDataSource(data=dict(value=[]))
inputs.on_change('data', recalculate_exptime)


eac_buttons = RadioButtonGroup(labels=EACS, active=0)
eac_buttons.js_on_event("button_click", CustomJS(args=dict(source=inputs, btn=eac_buttons), code="""
    source.data = { eacnum: [btn.active] }
"""))

newsnr  = Slider(title="Target SNR", value=10., start=0.1, end=100.0, step=0.1, ) 
newsnr_callback = CustomJS(args=dict(source=inputs), code="""
    source.data = { snr: [cb_obj.value] }
""")
newsnr.js_on_change("value", newsnr_callback)

star = Select(title="Template Star Spectrum", value="G2V star", 
                options=list(target_star.keys()), width=250) 
star_callback = CustomJS(args=dict(source=inputs), code="""
    source.data = { startemplate: [cb_obj.value] }
""")
star.js_on_change("value", star_callback)

planet = Select(title="Template Spectrum", value="Earth", 
                options=list(target_planet.keys()), width=250) 
planet_callback = CustomJS(args=dict(source=inputs), code="""
    source.data = { planettemplate: [cb_obj.value] }
""")
planet.js_on_change("value", planet_callback)

delta_mag  = Slider(title="delta Mag", value=15., start=10, end=30.0, step=0.1, ) 
dm_callback = CustomJS(args=dict(source=inputs), code="""
    source.data = { delta_mag: [cb_obj.value] }
""")
delta_mag.js_on_change("value", dm_callback)

texp_plot = figure(width=640, height=400, tools=("hover", "box_zoom", "wheel_zoom", "reset"), tooltips=[("@wave", "@flux")], toolbar_location="below")
texp_plot.line("wavelength", "exptime", source=obsdata)

controls = column(children=[eac_buttons,newsnr, star, planet, delta_mag], sizing_mode='fixed', max_width=300, width=300, height=700) 
#controls_tab = TabPanel(child=controls, title='Controls')
#plots_tab = TabPanel(child=texp_plot, title='Info')
l = layout([[controls, texp_plot]],sizing_mode='scale_width')

curdoc().add_root(l) 
curdoc().add_root(obsdata) 
