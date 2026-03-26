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
EACS = ["EAC1"]
#.         u'/U,  B,     g,     V      r,     Rc,    i,     Ic,    z',    y',    J,    H    
FILTERS = [0.365, 0.442, 0.475, 0.540, 0.617, 0.647, 0.752, 0.786, 0.866, 0.962, 1.22, 1.63]

class CoronImaging(pyEDITHETC):

    EACS = ["EAC1"]

    def __init__(self):
        self.intro = Div(text=f'<p>This Habworlds Coronagraphic Imaging ETC is powered by pyEDITH (E. Alei, M. Currie, C. Stark), v{pE.__version__}.</p><p>Selecting a planet will reset the default separation.</p>')

        super().__init__()

    def widget_setup(self):
        super().widget_setup()
        self.exp_plot = figure(height=480, title=f"", x_axis_label='microns', y_axis_label='Exposure Time (hr)', tools=("crosshair,pan,reset,save,box_zoom,wheel_zoom,hover"), tooltips=[("Wavelength (microns): ", "@wavelength"), ("Exposure Time (hr): ", "@exptime")], toolbar_location="below")
        self.exp_plot.scatter("wavelength", "exptime", source=self.obsdata)
        self.exp_panel = TabPanel(child=self.exp_plot, title='Exposure Time') #, width=800)


        self.snr_plot = figure(height=480, title=f"", x_axis_label='microns', y_axis_label='SNR', tools=("crosshair,pan,reset,save,box_zoom,wheel_zoom,hover"), tooltips=[("Wavelength (microns): ", "@wavelength"), ("SNR: ", "@snr")], toolbar_location="below")
        self.snr_plot.scatter("wavelength", "snr", source=self.obsdata)
        self.snr_panel = TabPanel(child=self.snr_plot, title='SNR') #, width=800)

    def load_initial(self):
        """
        Load initial parameters for the Imaging ETC
        """
        # observation parameters
        # set up wavelengths
        self.parameters["wavelength"] = 0.5
        self.parameters["snr"] = param_snr# the SNR you want for each spectral bin 
        self.parameters["CRb_multiplier"] = 2. # factor to multiply the background by (used for differential imaging)
        #self.parameters["photometric_aperture_radius"] = None#0.85 # radius of the photometric aperture in units of lambda/D
        self.parameters["psf_trunc_ratio"] = 0.3 # truncate the off-axis PSFs at this level 
        self.parameters["bandwidth"] = 0.2
        self.parameters["regrid_wavelength"] = False # set the flag to do this. We also need to specify a few other parameters.

        # The Astrophysical 
        # STAR
        #parameters["Lstar"] = 1. # luminosity of the star in solar luminosities
        self.parameters["distance"] = 10. # distance to the system in pc
        self.parameters["semimajor_axis"] = 1 # planetary separation in AU
        #parameters["Fp/Fs"] = FpFs # 1e-8 for testing (bright planet)
        # Note: we can work in either mag or flux units. Let's choose to work in flux units. 
        #Fstar_obs_10pc = compute_blackbody_photon_flux(5770, parameters["wavelength"] << u.um, parameters["distance"] << u.pc)

        self.parameters["stellar_radius"] = 1 # physical diameter of the star in Rsun units
        self.parameters["magV"] = 4.5 # ABMag
        # just this first time, because we need to load a planet to start with
        self.reflect_planet = self.target_planet["Earth"]["spectrum"]
        self.parameters["planetary_radius"] = self.target_planet["Earth"]["planetary_radius"]



        # PLANET
        # Most of these input spectra are albedo for a system. Actual conversion would require knowing the system separation and relative sizes of the star and planet. OR we make that an input parameter.
        self.parameters["delta_mag"] = 15# * u.ABmag
        # FpFs = target_planet["Earth"]["spectrum"](self.parameters["wavelength"] << u.um)
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

        # SCENE
        self.parameters["nzodis"] = 3. # number of zodis for exozodi estimate
        self.parameters["ra"] = 176.6292 # approximate ra of HD 102365. WARNING: do not use this number for science. 
        self.parameters["dec"] = -40.5003 # approximate dec of HD 102365. WARNING: do not use this number for science. 

        # Observatory parameters
        self.parameters["observing_mode"] = "IMAGER" # ETC should use IMAGER mode
        if "eacnum" in self.parameters:
            self.parameters["observatory_preset"] = EACS[self.parameters["eacnum"]]
        else:
            self.parameters["observatory_preset"] = "EAC1" # tells ETC to use EAC1 yaml files throughputs
        #self.parameters["npix_multiplier"] = np.ones_like(self.parameters["wavelength"]) # number of detector pixels per spectral bin
        self.parameters["noisefloor_PPF"] = 30 # post processing factor of 30 is a good realistic value for this

        # We need to run this, but parse_parameters also effectively strips out anything
        # pyEDITH itself doesn't use, and we're using the same dictionary for a lot of our own
        # parameters. 
        self.parameters = pE.parse_input.parse_parameters(self.parameters)

        # star must be loaded first, because planet flux is relative to star.
        self.load_star("G2V star")
        self.load_planet("Earth")
        self.recompute_star_flux()
        self.recompute_planet_flux()

        # this piece, alone, has to be created WITH some configured parameters.
        self.observatory_config = pE.parse_input.get_observatory_config(self.parameters)

        self.observatory = pE.ObservatoryBuilder.create_observatory(self.observatory_config)

        self.recalculate_exptime(ColumnDataSource(data={"scene": [True], "observatory": [True], "observation": [True]}))



    def do_recalculate_exptime(self, newvalues):
        global obsdata
        global exptime_compute

        wave_filters = []
        exptime_filters = []
        fpfs_filters = []
        obs_filters = []
        noise_hi = []
        noise_lo = []
        snr_filters = []
        self.update_calculation(newvalues)

        for filter in FILTERS:
            self.parameters["wavelength"] = [filter]
            self.update_calculation(ColumnDataSource(data={"scene": [True], "observatory": [True], "observation": [True]}))

            try:
                pE.calculate_exposure_time_or_snr(self.observation, self.scene, self.observatory, mode="exposure_time")
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
            wave_filters.append(self.observation.wavelength.to_value(u.um))
            exptime_filters.append(self.observation.exptime[0].to_value(u.hr))
            fpfs_filters.append(self.scene.Fp_over_Fs[0].value)
            obs_filters.append(obs[0].value)
            noise_hi.append(obs[0].value + noise[0].value/2.)
            noise_lo.append(obs[0].value - noise[0].value/2.)
            snr_filters.append(self.newsnr.value)
        exptime_filters = np.asarray(exptime_filters) << u.hr
        wave_filters = np.asarray(wave_filters) << u.um
        fpfs_filters = np.asarray(fpfs_filters) << u.dimensionless_unscaled
        obs_filters = np.asarray(obs_filters) << u.dimensionless_unscaled
        noise_hi = np.asarray(noise_hi) << u.dimensionless_unscaled
        noise_lo = np.asarray(noise_lo) << u.dimensionless_unscaled
        snr_filters = np.asarray(snr_filters) << u.dimensionless_unscaled
        print(exptime_filters)

        good = np.where(exptime_filters < 1e9 * u.s)[0] # there's no way we're doing anything that takes 100,000,000 seconds (3.169 years)
        #print(obs_filters, noise_hi)

        self.obsdata.data={"wavelength": wave_filters[good], "exptime": exptime_filters[good], "FpFs": fpfs_filters[good], 
                    "obs": obs_filters[good], "noise_hi": noise_hi[good], "noise_lo": noise_lo[good], "snr": snr_filters[good]}
    #print("New Data", obsdata.data)
        title_text = f"{self.planet.value} - {self.star.value} - {np.round(self.distance.value, decimals=2)} pc - {np.round(self.semimajor.value, decimals=2)} AU - SNR={np.round(self.newsnr.value, decimals=2)} - {EACS[self.eac_buttons.active]}"
        self.exp_plot.title.text =  title_text
        self.spec_plot.title.text = title_text

        self.exptime_compute.label = "Compute"



    def do_recalculate_snr(self, newvalues):
        self.update_calculation(newvalues)

        self.observation.obstime = (self.newexp.value * u.hr).to(u.s)

        wave_filters = []
        exptime_filters = []
        fpfs_filters = []
        obs_filters = []
        noise_hi = []
        noise_lo = []
        snr_filters = []
        self.update_calculation(newvalues)

        for filter in FILTERS:
            self.parameters["wavelength"] = [filter]
            self.update_calculation(ColumnDataSource(data={"scene": [True], "observatory": [True], "observation": [True]}))

            try:
                pE.calculate_exposure_time_or_snr(self.observation, self.scene, self.observatory, mode="signal_to_noise")
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
            obs, noise = pE.utils.synthesize_observation(self.observation.fullsnr,
                                                    self.scene, 
                                                    random_seed=None, # seed defaults to None
                                                    set_below_zero=0., # if the fake data falls below zero, set the data point as this. default = NaN
                                                    )

            print("Obs", obs)
            print("Noise", noise)
            wave_filters.append(self.observation.wavelength.to_value(u.um))
            exptime_filters.append(self.newexp.value)
            fpfs_filters.append(self.scene.Fp_over_Fs[0].value)
            obs_filters.append(obs[0].value)
            noise_hi.append(obs[0].value + noise[0].value/2.)
            noise_lo.append(obs[0].value - noise[0].value/2.)
            snr_filters.append(self.observation.fullsnr[0])
        exptime_filters = np.asarray(exptime_filters) << u.hr
        wave_filters = np.asarray(wave_filters) << u.um
        fpfs_filters = np.asarray(fpfs_filters) << u.dimensionless_unscaled
        obs_filters = np.asarray(obs_filters) << u.dimensionless_unscaled
        noise_hi = np.asarray(noise_hi) << u.dimensionless_unscaled
        noise_lo = np.asarray(noise_lo) << u.dimensionless_unscaled
        snr_filters = np.asarray(snr_filters) << u.dimensionless_unscaled
        print(exptime_filters)

        good = np.where(exptime_filters < 1e9 * u.s)[0] # there's no way we're doing anything that takes 100,000,000 seconds (3.169 years)
        #print(obs_filters, noise_hi)

        self.obsdata.data={"wavelength": wave_filters[good], "exptime": exptime_filters[good], "FpFs": fpfs_filters[good], 
                    "obs": obs_filters[good], "noise_hi": noise_hi[good], "noise_lo": noise_lo[good], "snr": snr_filters[good]}
    #print("New Data", obsdata.data)

        self.snr_plot.title.text =  f"{self.planet.value} - {self.star.value} - {np.round(self.distance.value, decimals=2)} pc - {np.round(self.semimajor.value, decimals=2)} AU - Exptime={np.round(self.newexp.value, decimals=2)} hrs - {EACS[self.eac_buttons.active]}"
        self.spec_plot.title.text =  f"{self.planet.value} - {self.star.value} - {np.round(self.distance.value, decimals=2)} pc - {np.round(self.semimajor.value, decimals=2)} AU - Exptime={np.round(self.newexp.value, decimals=2)} hrs - {EACS[self.eac_buttons.active]}"

        self.snr_compute.label = "Compute"






coron_imaging = CoronImaging()


coron_imaging.load_initial()

