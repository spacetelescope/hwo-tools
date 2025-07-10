import os
from functools import partial

import numpy as np
import scipy as sc
from bokeh.plotting import figure
from bokeh.embed import components
from bokeh.models import ColumnDataSource, Paragraph, Range1d, RadioButtonGroup, SetValue
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Slider, Div, Select, Button
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.layouts import row, column, layout
from bokeh.io import curdoc, output_file

import astropy.units as u
import astropy.constants as c
import synphot as syn
import stsynphot as stsyn

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
2. I don't think the plot data is actually reloading each time it runs.
3. Implement a "Calculate" button because PyEdith takes time, instead of any change to data.
4. Implement the load_star and load_planet features
5. Get the hover-text working
'''

param_snr=10
EACS = ["EAC1", "EAC2", "EAC3"]

target_planet, target_star = catalog.load_catalog()

parameters = {}
scene = pE.AstrophysicalScene()
observation = pE.Observation() # define the observation object
observatory = None # this piece, alone, has to be created WITH some configured parameters. So that's done in load_initial()
obsdata = ColumnDataSource(data=dict(wavelength=[], exptime=[]))
inputs = ColumnDataSource(data=dict())

texp_plot = figure(width=640, height=400, title=f"", x_axis_label=r'microns $$\{mu}$$ m', y_axis_label='Exposure Time (hr)', tools=("hover", "box_zoom", "wheel_zoom", "reset"), tooltips=[("@wavelength", "@exptime")], toolbar_location="below")
texp_plot.line("wavelength", "exptime", source=obsdata)

compute = Button(label="Calculate", button_type="primary")
# Can't set up the callback here because we need to define its callback (recalculate_exptime) first.

def compute_blackbody_photon_flux(temp, wavelengths, dist):
    """Generate photon flux density (photon/s/cm^2/um) for a blackbody at 1 cm^2."""
    bb = SourceSpectrum(BlackBodyNorm1D, temperature=temp)
    

    flux_photlam = syn.units.convert_flux(wavelengths, bb(wavelengths), u.photon / (u.s * u.cm**2 * u.nm)) #.value * 1/(u.s * u.cm**2 * u.AA) # photons/s/cm2/A
    #flux_photlam = flux_photlam.to(1 / (u.s * u.cm**2 * u.nm)) # convert to 1/nm
    flux_photlam = flux_photlam * (1000*u.pc/dist)**2

    return flux_photlam

def update_scene():
    global parameters
    global scene
    #parameters.update(updates)

    scene.load_configuration(parameters)
    scene.calculate_zodi_exozodi(parameters)
    scene.validate_configuration()

    return scene

def update_observation():
    global parameters
    global observation
    #parameters.update(updates)

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


def load_initial():
    """
    Load initial parameters for the ETC on load
    """
    global parameters
    global observatory
    global reflect_planet

    # observation parameters
    # set up wavelengths
    parameters["wavelength"] = np.linspace(0.35, 1.7, 1000)
    parameters["nlambd"] = len(parameters["wavelength"]) # number of wavelengths
    parameters["snr"] = param_snr * np.ones_like(parameters["wavelength"]) # the SNR you want for each spectral bin 
    parameters["CRb_multiplier"] = 2. # factor to multiply the background by (used for differential imaging)
    #parameters["photometric_aperture_radius"] = None#0.85 # radius of the photometric aperture in units of lambda/D
    parameters["psf_trunc_ratio"] = 0.3 # truncate the off-axis PSFs at this level 

    parameters["regrid_wavelength"] = True # set the flag to do this. We also need to specify a few other parameters.
    parameters["spectral_resolution"] = np.array([70, 140, 90]) #np.array([140])  # we're going to define three spectral channels. These are the spectral resolutions for each channel. i.e. all spectral bins in a given channel will have a fixed resolution.
    parameters["channel_bounds"] = np.array([0.5, 1.0]) #np.array([]) # specify the boundaries between the channels in um


    # The Astrophysical 
    # STAR
    parameters["Lstar"] = 1. # luminosity of the star in solar luminosities
    parameters["distance"] = 10. # distance to the system in pc
    parameters["separation"] = 0.1 # planetary separation in arcsec
    #parameters["Fp/Fs"] = FpFs # 1e-8 for testing (bright planet)
    parameters["Fp_min/Fs"] = 1e-11 # minimum detectable planet 
    # Note: we can work in either mag or flux units. Let's choose to work in flux units. 
    Fstar_obs_10pc = compute_blackbody_photon_flux(5770, parameters["wavelength"] << u.um, parameters["distance"] << u.pc)

    parameters["stellar_angular_diameter"] = 1.5e-4 # angular diameter of the star = 0.01 lam/D == about 1.5e-4 arcsec for lam=0.5 um and D~6.5 m
    parameters["magnitude"] = 4.5
    # just this first time, because we need to load a planet to start with
    reflect_planet = target_planet["Earth"]["spectrum"]
    parameters["planetary_radius"] = target_planet["Earth"]["planetary_radius"]

    load_star("G2V star")

    #parameters["Fstar_10pc"] = Fstar_obs_10pc.value # Fstar modeled as a blackbody (see above)
    #parameters["FstarV_10pc"] = v_band_flux.value # V-band flux of the star; you don't have to pre-calculate this, ETC will do it for you
    #print("Star", parameters["Fstar_10pc"])



    # PLANET
    # Most of these input spectra are albedo for a system. Actual conversion would require knowing the system separation and relative sizes of the star and planet. OR we make that an input parameter.
    parameters["delta_mag"] = 15 * u.ABmag
    # FpFs = target_planet["Earth"]["spectrum"](parameters["wavelength"] << u.um)
    # # normalize
    # contrast = 1e-8
    # FpFs = FpFs/np.max(FpFs)
    # Fplan_obs = Fstar_obs_10pc*contrast
    # FpFs = (Fplan_obs / Fstar_obs_10pc).value

    # from matplotlib import pyplot as plt
    # plt.plot(parameters["wavelength"] * u.um, Fstar_obs_10pc)
    # plt.plot(parameters["wavelength"] * u.um, Fplan_obs)
    # plt.yscale("log")

    # plt.show()

    #FpFs = FpFs * np.median(Fstar_obs_10pc) - syn.units.convert_flux(parameters["wavelength"], FpFs, parameters["delta_mag"]
    #print("Planet", FpFs, Fplan_obs)
    #flux_planet = Fstar_obs_10pc * FpFs * parameters["delta_mag"]
    #parameters["F0"] = flux_planet.value

    load_planet("Earth")
    #print("Star and Planet", parameters["Fstar_10pc"], parameters["Fp/Fs"])

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

    # this piece, alone, has to be created WITH some configured parameters.
    observatory_config = pE.parse_input.get_observatory_config(parameters)

    observatory = pE.ObservatoryBuilder.create_observatory(observatory_config)


    recalculate_exptime(ColumnDataSource(data={"scene": [True], "observatory": [True], "observation": [True]}))

def recompute_planet_flux():
    global parameters
    global reflect_planet
    solid_angle = parameters["planetary_radius"]**2/(4 * (parameters["separation"]*parameters["distance"]*1.5e8)**2) #Momentarily put both in km. pi cancels out of top and bottom. 
    flux_planet = parameters["Fstar_10pc"] * solid_angle * reflect_planet(parameters["wavelength"] << u.micron)
    parameters["F0"] = flux_planet.value
    parameters["Fp/Fs"] = (flux_planet / parameters["Fstar_10pc"]).value
    separation.value = parameters["separation"] # Make sure it matches what was used


    print("Star:", parameters["Fstar_10pc"])
    print("Planet:", reflect_planet(parameters["wavelength"] << u.micron), solid_angle)
    print(parameters["separation"], parameters["distance"])

def load_planet(sourceID):
    global parameters
    global reflect_planet
    # this is an albedo; the amount of incident flux received at that distance
    reflect_planet = target_planet[sourceID]["spectrum"]
    parameters["planetary_radius"] = target_planet[sourceID]["planetary_radius"]
    parameters["separation"] = .1 * target_planet[sourceID]["separation"]
    separation.value = .1 * target_planet[sourceID]["separation"]


    recompute_planet_flux() # trigger a recomputation of the planetary flux

def load_star(sourceID):
    global parameters
    star = target_star[sourceID]["spectrum"]
    parameters["magnitude"] = target_star[sourceID]["magnitude"] # Johnson V magnitude, specifically.
    magnitude.value = target_star[sourceID]["magnitude"]
    parameters["current_star"] = star

    recompute_star_flux()

def recompute_star_flux():

    bp = stsyn.band("johnson,v")
    #print(bp)
    #print(parameters["current_star"])

    new_star = parameters["current_star"].normalize(parameters["magnitude"] * u.ABmag, band=bp)
    flux = new_star(parameters["wavelength"]<< u.micron)
    magnitude.value = parameters["magnitude"] # make sure it matches what was used

    parameters["Fstar_10pc"] = syn.units.convert_flux(parameters["wavelength"], flux, u.photon / (u.s * u.cm**2 * u.nm)).value

    recompute_planet_flux() # trigger a recomputation of the planetary flux

    #print("Star Flux", flux)

# Ordinarily, these would be separate, but at the moment all changes here would seem to affect observatory, observation, and scene
# it is, particularly, unclear what 
def update_calculation(newvalues):
    global parameters
    global observation
    global scene
    global observatory
    print("------------------------------------")
    print(newvalues.data)

    if "new_star" in newvalues.data:
        print("Changed star")
        load_star(newvalues.data["new_star"][0])
        del newvalues.data["new_star"] # consume the new value
    if "new_magnitude" in newvalues.data:
        print("Changed stellar magnitude")
        parameters["magnitude"] = newvalues.data["new_magnitude"][0]
        recompute_star_flux()
        del newvalues.data["new_magnitude"] # consume the new value
    if "new_planet" in newvalues.data:
        print("Changed planet")
        load_planet(newvalues.data["new_planet"][0])
        del newvalues.data["new_planet"]
    if "new_separation" in newvalues.data:
        print("Changed Separation")
        parameters["separation"] = newvalues.data["new_separation"][0]
        recompute_planet_flux()
        del newvalues.data["new_separation"]
    if "new_snr" in newvalues.data:
        print("Changed SNR")
        parameters["snr"] = newvalues.data["new_snr"][0] * np.ones_like(parameters["wavelength"])
        del newvalues.data["new_snr"]
    if "new_eac" in newvalues.data:
        parameters["observatory_preset"] = newvalues.data["new_eac"][0]
        print("Changed EAC")
        del newvalues.data["new_eac"]
    if "new_diameter" in newvalues.data:
        print("Changed Diameter")
        parameters["diameter"] = newvalues.data["new_diameter"][0]
        del newvalues.data["new_diameter"]
    else:
        parameters["observatory_preset"] = "EAC1" # tells ETC to use EAC1 yaml files throughputs

    if "observation" in newvalues.data and newvalues.data["observation"][0]:
        print("Rerun observation...")
        observation = update_observation()

    if "scene" in newvalues.data and newvalues.data["scene"][0]:
        print("Rerun scene")
        scene = update_scene()
        if parameters["regrid_wavelength"] is True:
            scene.regrid_spectra(parameters, observation)

    pE.ObservatoryBuilder.configure_observatory(
        observatory, parameters, observation, scene
    )
    observatory.validate_configuration()
    #print(observatory.telescope.__dict__)
    #print_observatory(observatory)


    return observatory, scene, observation

def do_recalculate_exptime(newvalues):

    observatory, scene, observation = update_calculation(newvalues)

    pE.calculate_exposure_time_or_snr(observation, scene, observatory, verbose=True)

    good = np.where(observation.exptime < 1e8 * u.s)

    obsdata.data={"wavelength": observation.wavelength[good], "exptime": observation.exptime[good].to(u.hr)}
    #print("New Data", obsdata.data)
    texp_plot.title.text =  f"{planet.value} - {star.value} - SNR={np.round(newsnr.value, decimals=2)} - {EACS[eac_buttons.active]}"

    compute.label = "Compute"


def recalculate_exptime(newvalues):
    """
    The trick here is that Bokeh only synchronizes the calls at the end of a function, so
    if I want to change the button to "Please wait..." and THEN have it calculate, I have
    to make this call, and have this call fire off another callback with add_next_tick_callback.

    Parameters
    ----------
    newvalues : _type_
        _description_
    """
    # so we can redo all data
    global obsdata
    global compute

    compute.label = "Please Wait..."
    curdoc().add_next_tick_callback(partial(do_recalculate_exptime, newvalues))


compute.on_click(partial(recalculate_exptime, inputs))

def recalculate_snr(newvalues):
    # so we can redo all data
    global obsdata
    observatory, scene, observation = update_calculation(newvalues)

    pE.calculate_exposure_time_or_snr(observation, scene, observatory, mode="signal_to_ noise", verbose=False)

    good = np.where(~np.isinf(observation.fullsnr))

    obsdata.data={"wavelength": observation.wavelength[good], "exptime": observation.fullsnr[good]}

eac_buttons = RadioButtonGroup(labels=EACS, active=0)

def eac_callback(attr, old, new):
    global inputs
    print(attr, old, new)
    inputs.data.update({"new_eac": [EACS[new]], "observatory": [True]})
eac_buttons.on_change("active", eac_callback)

newsnr  = Slider(title="Target SNR", value=10., start=0.1, end=100.0, step=0.1, ) 
def snr_callback(attr, old, new):
    global inputs
    print(attr, old, new)
    inputs.data.update({"new_snr": [new], "observation": [True]})
newsnr.on_change("value", snr_callback)

newdiameter  = Slider(title="Mirror Diameter", value=7., start=5, end=15, step=0.1, ) 
def diameter_callback(attr, old, new):
    global inputs
    print(attr, old, new)
    inputs.data.update({"new_diameter": [new], "observatory": [True]})
newdiameter.on_change("value", diameter_callback)

star = Select(title="Template Star Spectrum", value="G2V star", 
                options=list(target_star.keys()), width=250) 
def star_callback(attr, old, new):
    global inputs
    print(attr, old, new)
    inputs.data.update({"new_star": [new], "scene": [True]})
star.on_change("value", star_callback)

magnitude  = Slider(title="Stellar Magnitude", value=4.5, start=3, end=20.0, step=0.1) 
def magnitude_callback(attr, old, new):
    global inputs
    print(attr, old, new)
    inputs.data.update({"new_magnitude": [new], "scene": [True]})
magnitude.on_change("value", magnitude_callback)

planet = Select(title="Template Planet Spectrum", value="Earth", 
                options=list(target_planet.keys()), width=250)
def planet_callback(attr, old, new):
    global inputs
    print(attr, old, new)
    inputs.data.update({"new_planet": [new], "scene": [True]})
planet.on_change("value", planet_callback)

separation  = Slider(title="Separation (arcsec @ 10pc)", value=0.1, start=0.01, end=0.5, step=0.01, ) 
def separation_callback(attr, old, new):
    global inputs
    print(attr, old, new)
    inputs.data.update({"new_separation": [new], "scene": [True]})
separation.on_change("value", separation_callback)

delta_mag  = Slider(title="delta Mag", value=15., start=10, end=30.0, step=0.1, ) 
def dmag_callback(attr, old, new):
    global inputs
    print(attr, old, new)
    inputs.data.update({"new_dMag": [new], "scene": [True]})
delta_mag.on_change("value", dmag_callback)


load_initial()


controls = column(children=[eac_buttons, newdiameter, newsnr, star, magnitude, planet, separation, compute], sizing_mode='fixed', max_width=300, width=300, height=700) 
#controls_tab = TabPanel(child=controls, title='Controls')
#plots_tab = TabPanel(child=texp_plot, title='Info')
l = layout([[controls, texp_plot]],sizing_mode='scale_width')

curdoc().theme = 'dark_minimal'
curdoc().add_root(l) 
curdoc().add_root(obsdata) 
