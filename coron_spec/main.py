import os
from functools import partial

import base64
import datetime
import asdf
import yaml
import copy

import numpy as np
import scipy as sc
from bokeh.plotting import figure
from bokeh.embed import components
from bokeh.models import ColumnDataSource, Paragraph, Range1d, RadioGroup, RadioButtonGroup, SetValue
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Slider, Div, Select, Button, FileInput
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.layouts import row, column, layout
from bokeh.io import curdoc, output_file

import astropy.units as u
import astropy.constants as c
import synphot as syn
import stsynphot as stsyn

import pyEDITH as pE
import eacy
from synphot import SourceSpectrum, SpectralElement, Observation
from synphot.models import BlackBodyNorm1D, Empirical1D
from synphot import units as synphot_u

import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/..")

from common.pyedith_etc_common import pyEDITHETC

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


class CoronSpec(pyEDITHETC):
    
    EACS = ["EAC1"]

    def __init__(self):
        self.intro = Div(text=f'<p>This Habworlds Coronagraphic Spectroscopy ETC is powered by pyEDITH (E. Alei, M. Currie, C. Stark), v{pE.__version__}.</p><p>Selecting a planet will reset the default separation.</p>')
        super().__init__()

    def load_initial(self):
        """
        Load initial parameters for the Spectroscopic ETC
        """

        # observation parameters
        # set up wavelengths
        self.parameters["wavelength"] = np.linspace(0.35, 1.71, 1000)
        self.parameters["nlambd"] = len(self.parameters["wavelength"]) # number of wavelengths
        self.parameters["snr"] = param_snr * np.ones_like(self.parameters["wavelength"]) # the SNR you want for each spectral bin 
        self.parameters["CRb_multiplier"] = 2. # factor to multiply the background by (used for differential imaging)
        #self.parameters["photometric_aperture_radius"] = None#0.85 # radius of the photometric aperture in units of lambda/D
        self.parameters["psf_trunc_ratio"] = 0.3 # truncate the off-axis PSFs at this level 

        self.parameters["regrid_wavelength"] = True # set the flag to do this. We also need to specify a few other parameters.
        self.parameters["spectral_resolution"] = np.array([70, 140, 90]) #np.array([140])  # we're going to define three spectral channels. These are the spectral resolutions for each channel. i.e. all spectral bins in a given channel will have a fixed resolution.
        self.parameters["channel_bounds"] = np.array([0.5, 1.0]) #np.array([]) # specify the boundaries between the channels in um
        self.parameters["lam_low"] = [0.36, 0.5, 1.0]
        self.parameters["lam_high"] = [0.5, 1.0, 1.7]

        # The Astrophysical 
        # STAR
        #self.parameters["Lstar"] = 1. # luminosity of the star in solar luminosities
        self.parameters["distance"] = 10. # distance to the system in pc
        self.parameters["semimajor_axis"] = 1 # planetary separation in AU
        #self.parameters["Fp/Fs"] = FpFs # 1e-8 for testing (bright planet)


        #parameters["Fstar_10pc"] = Fstar_obs_10pc.value # Fstar modeled as a blackbody (see above)
        #parameters["FstarV_10pc"] = v_band_flux.value # V-band flux of the star; you don't have to pre-calculate this, ETC will do it for you
        #print("Star", parameters["Fstar_10pc"])



        # PLANET
        # Most of these input spectra are albedo for a system. Actual conversion would require knowing the system separation and relative sizes of the star and planet. OR we make that an input parameter.
        self.parameters["delta_mag"] = 15# * u.ABmag
        # FpFs = self.target_planet["Earth"]["spectrum"](parameters["wavelength"] << u.um)
        # # normalize
        # contrast = 1e-8
        # FpFs = FpFs/np.max(FpFs)
        # Fplan_obs = Fstar_obs_10pc*contrast
        # FpFs = (Fplan_obs / Fstar_obs_10pc).value

        # from matplotlib import pyplot as plt
        # plt.plot(self.parameters["wavelength"] * u.um, Fstar_obs_10pc)
        # plt.plot(self.parameters["wavelength"] * u.um, Fplan_obs)
        # plt.yscale("log")

        # plt.show()

        #FpFs = FpFs * np.median(Fstar_obs_10pc) - syn.units.convert_flux(self.parameters["wavelength"], FpFs, self.parameters["delta_mag"]
        #print("Planet", FpFs, Fplan_obs)
        #flux_planet = Fstar_obs_10pc * FpFs * self.parameters["delta_mag"]
        #self.parameters["F0"] = flux_planet.value

        self.load_star("G2V star")
        self.load_planet("Earth")
        self.recompute_star_flux()
        self.recompute_planet_flux()
        #print("Star and Planet", self.parameters["Fstar_10pc"], self.parameters["Fp/Fs"])

        # SCENE
        self.parameters["nzodis"] = 3. # number of zodis for exozodi estimate
        self.parameters["ra"] = 176.6292 # approximate ra of HD 102365. WARNING: do not use this number for science. 
        self.parameters["dec"] = -40.5003 # approximate dec of HD 102365. WARNING: do not use this number for science. 

        # Observatory parameters
        self.parameters["observing_mode"] = "IFS" # ETC should use IFS mode
        if "eacnum" in self.parameters:
            self.parameters["observatory_preset"] = self.EACS[self.parameters["eacnum"]]
        else:
            self.parameters["observatory_preset"] = "EAC1" # tells ETC to use EAC1 yaml files throughputs
        self.parameters["IFS_eff"]  = 1. # extra throughput of the IFS 
        self.parameters["npix_multiplier"] = np.ones_like(self.parameters["wavelength"]) # number of detector pixels per spectral bin
        self.parameters["noisefloor_PPF"] = 30 # post processing factor of 30 is a good realistic value for this

        # this piece, alone, has to be created WITH some configured parameters.
        self.observatory_config = pE.parse_input.get_observatory_config(self.parameters)

        self.observatory = pE.ObservatoryBuilder.create_observatory(self.observatory_config)

        self.recalculate_exptime(ColumnDataSource(data={"scene": [True], "observatory": [True], "observation": [True]}))

    def do_recalculate_exptime(self, newvalues):
        global obsdata
        global exptime_compute

        self.update_calculation(newvalues)

        try:
            pE.calculate_exposure_time_or_snr(self.observation, self.scene, self.observatory)
        except UnboundLocalError:
            self.warning.text = "<p style='color:Tomato;'>ERROR: Inputs out of bounds. Try again</p>"
            self.exptime_compute.label = "Compute"
            self.obsdata.data={"wavelength": [], "exptime": [], "FpFs": [], "obs": [], "noise_hi": [], "noise_lo": [], "snr": []}

            return
        #print("SNR", newsnr.value * np.ones_like(observation.wavelength.value))
        #print("Exptime", observation.exptime)
        if any(np.isinf(self.observation.exptime)):
            self.warning.text = "<p style='color:Gold;'>WARNING: Planet outside OWA or inside IWA. Hardcoded infinity results.</p>"
        else:
            self.warning.text = "<p></p>"
        obs, noise = pE.utils.synthesize_observation(self.newsnr.value * np.ones_like(self.observation.wavelength.value),
                                                self.scene, 
                                                random_seed=None, # seed defaults to None
                                                set_below_zero=0., # if the fake data falls below zero, set the data point as this. default = NaN
                                                )

        print("Obs", obs)
        print("Noise", noise)

        good = np.where(self.observation.exptime < 1e8 * u.s) # there's no way we're doing anything that takes 100,000,000 seconds (3.169 years)

        self.obsdata.data={"wavelength": self.observation.wavelength[good], "exptime": self.observation.exptime[good].to(u.hr), "FpFs": self.scene.Fp_over_Fs[good], 
                    "obs": obs[good], "noise_hi": obs[good] + noise[good]/2., "noise_lo": obs[good] - noise[good]/2., "snr": self.newsnr.value * np.ones_like(self.observation.wavelength[good].value)}
        #print("New Data", obsdata.data)
        title_text = f"{self.planet.value} - {self.star.value} - {np.round(self.distance.value, decimals=2)} pc - {np.round(self.semimajor.value, decimals=2)} AU - SNR={np.round(self.newsnr.value, decimals=2)} - {self.EACS[self.eac_buttons.active]}"
        self.exp_plot.title.text =  title_text
        self.spec_plot.title.text = title_text

        self.exptime_compute.label = "Compute"


    def do_recalculate_snr(self, newvalues):
        self.update_calculation(newvalues)

        self.observation.obstime = (self.newexp.value * u.hr).to(u.s)
        print((self.newexp.value * u.hr).to(u.s))

        try:
            pE.calculate_exposure_time_or_snr(self.observation, self.scene, self.observatory, mode="signal_to_noise")
        except UnboundLocalError:
            self.warning.text = "<p style='color:Tomato;'>ERROR: Inputs out of bounds. Try again</p>"
            self.snr_compute.label = "Compute"
            self.obsdata.data={"wavelength": [], "exptime": [], "FpFs": [], "obs": [], "noise_hi": [], "noise_lo": [], "snr": []}

            return
        #print("SNR", newsnr.value * np.ones_like(observation.wavelength.value))
        #print("Exptime", observation.exptime)
        if any(np.isinf(self.observation.exptime)):
            self.warning.text = "<p style='color:Gold;'>WARNING: Planet outside OWA or inside IWA. Hardcoded infinity results.</p>"
        else:
            self.warning.text = "<p></p>"
        obs, noise = pE.utils.synthesize_observation(self.observation.fullsnr,
                                                self.scene, 
                                                random_seed=None, # seed defaults to None
                                                set_below_zero=0., # if the fake data falls below zero, set the data point as this. default = NaN
                                                )

        print("Obs", obs)
        print("Noise", noise)

        good = np.where(self.observation.exptime < 1e8 * u.s) # there's no way we're doing anything that takes 100,000,000 seconds (3.169 years)

        self.obsdata.data={"wavelength": self.observation.wavelength, "exptime": self.newexp.value * np.ones_like(self.observation.wavelength.value), "FpFs": self.scene.Fp_over_Fs, "obs": obs, "noise_hi": obs + noise/2., "noise_lo": obs - noise/2., "snr": self.observation.fullsnr}
        #print("New Data", obsdata.data)
        title_text = f"{self.planet.value} - {self.star.value} - {np.round(self.distance.value, decimals=2)} pc - {np.round(self.semimajor.value, decimals=2)} AU - Exptime={np.round(self.newexp.value, decimals=2)} hrs - {self.EACS[self.eac_buttons.active]}"
        self.snr_plot.title.text =  title_text
        self.spec_plot.title.text = title_text

        self.snr_compute.label = "Compute"


coron_spec = CoronSpec()


coron_spec.load_initial()
