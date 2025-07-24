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
EACS = ["EAC1"]

target_planet, target_star = catalog.load_catalog()

parameters = {}
scene = pE.AstrophysicalScene()
observation = pE.Observation() # define the observation object
observatory = None # this piece, alone, has to be created WITH some configured parameters. So that's done in load_initial()
obsdata = ColumnDataSource(data=dict(wavelength=[], exptime=[], FpFs=[], obs=[], noise_hi=[], noise_lo=[], snr=[]))
inputs = ColumnDataSource(data=dict())

texp_plot = figure(height=480, title=f"", x_axis_label='microns', y_axis_label='Exposure Time (hr)', tools=("crosshair,pan,reset,save,box_zoom,wheel_zoom,hover"), tooltips=[("Wavelength (microns): ", "@wavelength"), ("Exposure Time (hr): ", "@exptime")], toolbar_location="below")
texp_plot.line("wavelength", "exptime", source=obsdata)

exp_panel = TabPanel(child=texp_plot, title='Exposure Time') #, width=800)

snr_plot = figure(height=480, title=f"", x_axis_label='microns', y_axis_label='Fp/Fs', tools=("crosshair,pan,reset,save,box_zoom,wheel_zoom,hover"), tooltips=[("Wavelength (microns): ", "@wavelength"), ("Fp/Fs: ", "@FpFs"), ("SNR: ", "@snr")], toolbar_location="below")
snr_plot.line("wavelength", "FpFs", source=obsdata)
snr_plot.scatter('wavelength', 'obs', source=obsdata, fill_color='#B4D9FF', line_color='black', size=8, name='snr_plot_circle_hover') 
snr_plot.segment('wavelength', 'noise_hi', 'wavelength', 'noise_lo', source=obsdata, line_width=1, line_color='#82AFF6', line_alpha=0.5)
snr_panel = TabPanel(child=snr_plot, title='Spectrum') #, width=800)

compute = Button(label="Calculate", button_type="primary")
# Can't set up the callback here because we need to define its callback (recalculate_exptime) first.

warning = Div(text='<p></p>')

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
    parameters["wavelength"] = np.linspace(0.35, 1.71, 1000)
    parameters["nlambd"] = len(parameters["wavelength"]) # number of wavelengths
    parameters["snr"] = param_snr * np.ones_like(parameters["wavelength"]) # the SNR you want for each spectral bin 
    parameters["CRb_multiplier"] = 2. # factor to multiply the background by (used for differential imaging)
    #parameters["photometric_aperture_radius"] = None#0.85 # radius of the photometric aperture in units of lambda/D
    parameters["psf_trunc_ratio"] = 0.3 # truncate the off-axis PSFs at this level 

    parameters["regrid_wavelength"] = True # set the flag to do this. We also need to specify a few other parameters.
    parameters["spectral_resolution"] = np.array([70, 140, 90]) #np.array([140])  # we're going to define three spectral channels. These are the spectral resolutions for each channel. i.e. all spectral bins in a given channel will have a fixed resolution.
    parameters["channel_bounds"] = np.array([0.5, 1.0]) #np.array([]) # specify the boundaries between the channels in um
    parameters["lam_low"] = [0.36, 0.5, 1.0]
    parameters["lam_high"] = [0.5, 1.0, 1.7]

    # The Astrophysical 
    # STAR
    #parameters["Lstar"] = 1. # luminosity of the star in solar luminosities
    parameters["distance"] = 10. # distance to the system in pc
    parameters["semimajor_axis"] = 1 # planetary separation in AU
    #parameters["Fp/Fs"] = FpFs # 1e-8 for testing (bright planet)
    # Note: we can work in either mag or flux units. Let's choose to work in flux units. 
    Fstar_obs_10pc = compute_blackbody_photon_flux(5770, parameters["wavelength"] << u.um, parameters["distance"] << u.pc)

    parameters["stellar_radius"] = 1 # physical diameter of the star in Rsun units
    parameters["magV"] = 4.5 # ABMag
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
    solid_angle = parameters["planetary_radius"]**2/(4 * (parameters["semimajor_axis"]*1.5e8)**2) #Momentarily put both in km. pi cancels out of top and bottom. 
    flux_planet = parameters["FstarV_10pc"] * solid_angle * reflect_planet(parameters["wavelength"] << u.micron)
    parameters["F0"] = flux_planet.value
    parameters["Fp/Fs"] = (flux_planet / parameters["FstarV_10pc"]).value
    semimajor.value = parameters["semimajor_axis"] # Make sure it matches what was used


    print("Star:", parameters["FstarV_10pc"])
    print("Planet:", reflect_planet(parameters["wavelength"] << u.micron), solid_angle)
    print(parameters["semimajor_axis"], parameters["distance"])

def load_planet(sourceID):
    global parameters
    global reflect_planet
    # this is an albedo; the amount of incident flux received at that distance
    reflect_planet = target_planet[sourceID]["spectrum"]
    parameters["planetary_radius"] = target_planet[sourceID]["planetary_radius"]
    parameters["semimajor_axis"] = target_planet[sourceID]["semimajor_axis"]
    semimajor.value = target_planet[sourceID]["semimajor_axis"]


    recompute_planet_flux() # trigger a recomputation of the planetary flux

def load_star(sourceID):
    global parameters
    star = target_star[sourceID]["spectrum"]
    parameters["magV"] = target_star[sourceID]["magV"] # Johnson V magnitude, specifically.
    parameters["stellar_radius"] = target_star[sourceID]["stellar_radius"]
    #magnitude.value = target_star[sourceID]["magV"]
    parameters["current_star"] = star

    recompute_star_flux()

def recompute_star_flux():

    bp = stsyn.band("johnson,v")
    #print(bp)
    #print(parameters["current_star"])

    # we do not move the star from 10 pc, we merely provide the magnitude (and flux) at 10 pc and pyEDITH does the rest
    new_star = parameters["current_star"].normalize(parameters["magV"] * u.ABmag, band=bp)
    flux = new_star(parameters["wavelength"]<< u.micron)
    #magnitude.value = parameters["magV"] # make sure it matches what was used

    parameters["Fstar_10pc"] = syn.units.convert_flux(parameters["wavelength"], flux, u.photon / (u.s * u.cm**2 * u.nm)).value
    # get the flux at 10 pc in the V band
    parameters["FstarV_10pc"] = syn.units.convert_flux(parameters["wavelength"], syn.Observation(new_star, bp).effstim(), u.photon / (u.s * u.cm**2 * u.nm)).value

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
        parameters["magV"] = newvalues.data["new_magnitude"][0]
        recompute_star_flux()
        del newvalues.data["new_magnitude"] # consume the new value
    if "new_planet" in newvalues.data:
        print("Changed planet")
        load_planet(newvalues.data["new_planet"][0])
        del newvalues.data["new_planet"]
    if "new_semimajor" in newvalues.data:
        print("Changed Semimajor Axis")
        parameters["semimajor_axis"] = newvalues.data["new_semimajor"][0]
        recompute_planet_flux()
        del newvalues.data["new_semimajor"]
    if "new_distance" in newvalues.data:
        print("Changed System distance")
        parameters["distance"] = newvalues.data["new_distance"][0]
        recompute_planet_flux()
        del newvalues.data["new_distance"]
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
        newvalues.data["observation"][0] = False

    if "scene" in newvalues.data and newvalues.data["scene"][0]:
        print("Rerun scene")
        scene = update_scene()
        if parameters["regrid_wavelength"] is True:
            scene.regrid_spectra(parameters, observation)
        newvalues.data["scene"][0] = False

    pE.ObservatoryBuilder.configure_observatory(
        observatory, parameters, observation, scene
    )
    observatory.validate_configuration()
    #print(observatory.telescope.__dict__)
    #print_observatory(observatory)


    return observatory, scene, observation

def do_recalculate_exptime(newvalues):

    observatory, scene, observation = update_calculation(newvalues)

    try:
        pE.calculate_exposure_time_or_snr(observation, scene, observatory, verbose=True)
    except UnboundLocalError:
        warning.text = "<p style='color:Tomato;'>ERROR: Inputs out of bounds. Try again</p>"
        compute.label = "Compute"
        obsdata.data={"wavelength": [], "exptime": [], "FpFs": [], "obs": [], "noise_hi": [], "noise_lo": [], "snr": []}

        return
    #print("SNR", newsnr.value * np.ones_like(observation.wavelength.value))
    #print("Exptime", observation.exptime)
    if any(np.isinf(observation.exptime)):
        warning.text = "<p style='color:Gold;'>WARNING: Planet outside OWA or inside IWA. Hardcoded infinity results.</p>"
    else:
        warning.text = "<p></p>"
    obs, noise = pE.utils.synthesize_observation(newsnr.value * np.ones_like(observation.wavelength.value),
                                             observation.exptime,
                                             observation.wavelength,
                                             observation,
                                             scene, 
                                             random_seed=None, # seed defaults to None
                                             set_below_zero=0., # if the fake data falls below zero, set the data point as this. default = NaN
                                             )

    print("Obs", obs)
    print("Noise", noise)

    good = np.where(observation.exptime < 1e8 * u.s) # there's no way we're doing anything that takes 100,000,000 seconds (3.169 years)

    obsdata.data={"wavelength": observation.wavelength[good], "exptime": observation.exptime[good].to(u.hr), "FpFs": scene.Fp_over_Fs[good], "obs": obs[good], "noise_hi": obs[good] + noise[good]/2., "noise_lo": obs[good] - noise[good]/2., "snr": newsnr.value * np.ones_like(observation.wavelength[good].value)}
    #print("New Data", obsdata.data)
    texp_plot.title.text =  f"{planet.value} - {star.value} - {np.round(distance.value, decimals=2)} pc - {np.round(semimajor.value, decimals=2)} AU - SNR={np.round(newsnr.value, decimals=2)} - {EACS[eac_buttons.active]}"
    snr_plot.title.text =  f"{planet.value} - {star.value} - {np.round(distance.value, decimals=2)} pc - {np.round(semimajor.value, decimals=2)} AU - SNR={np.round(newsnr.value, decimals=2)} - {EACS[eac_buttons.active]}"

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

    pE.calculate_exposure_time_or_snr(observation, scene, observatory, mode="signal_to_noise", verbose=False)

    good = np.where(~np.isinf(observation.fullsnr))

    obsdata.data={"wavelength": observation.wavelength[good], "exptime": observation.fullsnr[good]}


intro = Div(text='<p>This Habworlds Coronagraphic ETC is powered by PyEDITH (E. Alei, M. Currie, C. Stark).</p><p>Selecting a planet will reset the default separation.</p>')

info_panel = Div(sizing_mode="inherit", text="pyEDITH is a Python-based coronagraphic exposure time calculator built for the Habitable Worlds Observatory (HWO).<p>It is designed to simulate wavelength-dependent exposure times and SNR for both photometric and spectroscopic direct imaging observations. pyEDITH interfaces with engineering specifications defined by the HWO exploratory analytic cases, and allows the user to provide target system information, as well as alter observatory parameters for trade studies, to calculate synthetic HWO observations of Earth-like exoplanets. pyEDITH has heritage from the exposure time calculator built for the Altruistic Yield Optimizer (<a href='https://ui.adsabs.harvard.edu/abs/2014ApJ...795..122S/abstract'>C.C. Stark et al., 2014</a>), and has been validated against the AYO, exoSIMS, and EBS exposure time calculators.")
observation_tab = TabPanel(child=texp_plot, title='Observation') # , width=400)

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

# magnitude  = Slider(title="Stellar Magnitude (Johnson V)", value=4.5, start=0.0, end=20.0, step=0.1) 
# def magnitude_callback(attr, old, new):
#     global inputs
#     print(attr, old, new)
#     inputs.data.update({"new_magnitude": [new], "scene": [True]})
# magnitude.on_change("value", magnitude_callback)

distance  = Slider(title="Distance to System (pc)", value=10, start=1.4, end=100.0, step=0.1) 
def distance_callback(attr, old, new):
    global inputs
    print(attr, old, new)
    inputs.data.update({"new_distance": [new], "scene": [True]})
distance.on_change("value", distance_callback)

planet = Select(title="Template Planet Spectrum", value="Earth", 
                options=list(target_planet.keys()), width=250)
def planet_callback(attr, old, new):
    global inputs
    print(attr, old, new)
    inputs.data.update({"new_planet": [new], "scene": [True]})
planet.on_change("value", planet_callback)

semimajor = Slider(title="Semimajor Axis (AU)", value=0.1, start=0.01, end=10, step=0.01, ) 
def semimajor_callback(attr, old, new):
    global inputs
    print(attr, old, new)
    inputs.data.update({"new_semimajor": [new], "scene": [True]})
semimajor.on_change("value", semimajor_callback)

delta_mag = Slider(title="delta Mag", value=15., start=10, end=30.0, step=0.1, ) 
def dmag_callback(attr, old, new):
    global inputs
    print(attr, old, new)
    inputs.data.update({"new_dMag": [new], "scene": [True]})
delta_mag.on_change("value", dmag_callback)

info_panel = TabPanel(child=info_panel, title='Info') #, width=800)
load_initial()


controls = column(children=[intro, newdiameter, newsnr, star, distance, planet, semimajor, compute, warning], sizing_mode='fixed', width=320, height=480) 
#controls_tab = TabPanel(child=controls, title='Controls')
#plots_tab = TabPanel(child=texp_plot, title='Info')
outputs = Tabs(tabs=[ snr_panel, exp_panel, info_panel], sizing_mode="inherit")
plots = column(children=[outputs], sizing_mode='fixed', width=640, height=480)
l = layout([[controls, plots]],sizing_mode='fixed', width=960, height=480)

curdoc().theme = 'dark_minimal'
curdoc().add_root(l) 
curdoc().add_root(obsdata) 
