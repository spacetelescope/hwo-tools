import os
from functools import partial

import base64
import datetime

import numpy as np
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Paragraph, Range1d, RadioGroup, RadioButtonGroup, SetValue
from bokeh.models.widgets import Slider, Div, Select, Button, FileInput
from bokeh.models.layouts import TabPanel, Tabs
from bokeh.layouts import row, column, layout
from bokeh.io import curdoc

import astropy.units as u
import astropy.constants as c
import synphot as syn
import stsynphot as stsyn

import pyEDITH as pE

import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/..")

from common import catalog, pyedith_etc_common

param_snr=10

class CoronSpec(pyedith_etc_common.pyEDITHETC):
    # classmethods
    target_planet, target_star = catalog.load_catalog()
    EACS = ["EAC1"]

    def __init__(self):
        self.parameters = {}
        self.scene = pE.AstrophysicalScene()
        self.observation = pE.Observation() # define the observation object
        self.observatory = None # this piece, alone, has to be created WITH some configured parameters. So that's done in load_initial()
        self.obsdata = ColumnDataSource(data=dict(wavelength=[], exptime=[], FpFs=[], obs=[], noise_hi=[], noise_lo=[], snr=[]))
        self.inputs = ColumnDataSource(data=dict())

        self.widget_setup()
        self.tab_setup()
    
    def widget_setup(self):
        self.intro = Div(text=f'<p>This Habworlds Coronagraphic Spectroscopy ETC is powered by pyEDITH (E. Alei, M. Currie, C. Stark), v{pE.__version__}.</p><p>Selecting a planet will reset the default separation.</p>')

        self.exp_plot = figure(height=480, title=f"", x_axis_label='microns', y_axis_label='Exposure Time (hr)', tools=("crosshair,pan,reset,save,box_zoom,wheel_zoom,hover"), tooltips=[("Wavelength (microns): ", "@wavelength"), ("Exposure Time (hr): ", "@exptime")], toolbar_location="below")
        self.exp_plot.line("wavelength", "exptime", source=self.obsdata)
        self.exp_panel = TabPanel(child=self.exp_plot, title='Exposure Time') #, width=800)

        self.snr_plot = figure(height=480, title=f"", x_axis_label='microns', y_axis_label='SNR', tools=("crosshair,pan,reset,save,box_zoom,wheel_zoom,hover"), tooltips=[("Wavelength (microns): ", "@wavelength"), ("SNR: ", "@snr")], toolbar_location="below")
        self.snr_plot.line("wavelength", "snr", source=self.obsdata)
        self.snr_panel = TabPanel(child=self.snr_plot, title='SNR') #, width=800)


        self.spec_plot = figure(height=480, title=f"", x_axis_label='microns', y_axis_label='Fp/Fs', tools=("crosshair,pan,reset,save,box_zoom,wheel_zoom,hover"), tooltips=[("Wavelength (microns): ", "@wavelength"), ("Fp/Fs: ", "@FpFs"), ("SNR: ", "@snr")], toolbar_location="below")
        self.spec_plot.line("wavelength", "FpFs", source=self.obsdata)
        self.spec_plot.scatter('wavelength', 'obs', source=self.obsdata, fill_color='#B4D9FF', line_color='black', size=8, name='snr_plot_circle_hover') 
        self.spec_plot.segment('wavelength', 'noise_hi', 'wavelength', 'noise_lo', source=self.obsdata, line_width=1, line_color='#82AFF6', line_alpha=0.5)
        self.spec_panel = TabPanel(child=self.spec_plot, title='Spectrum') #, width=800)

        self.exptime_compute = Button(label="Calculate", button_type="primary")
        # Can't set up the callback here because we need to define its callback (recalculate_exptime) first.
        self.exptime_compute.on_click(partial(self.recalculate_exptime, self.inputs))

        self.snr_compute = Button(label="Calculate", button_type="primary")
        # Can't set up the callback here because we need to define its callback (recalculate_snr) first.
        self.snr_compute.on_click(partial(self.recalculate_snr, self.inputs))

        self.warning = Div(text='<p></p>')

        self.info_panel = Div(sizing_mode="inherit", text="pyEDITH is a Python-based coronagraphic exposure time calculator built for the Habitable Worlds Observatory (HWO)." +
                    "<p>It is designed to simulate wavelength-dependent exposure times and SNR for both photometric and spectroscopic direct imaging observations. pyEDITH interfaces with engineering specifications defined by the HWO exploratory analytic cases, and allows the user to provide target system information, as well as alter observatory parameters for trade studies, to calculate synthetic HWO observations of Earth-like exoplanets. pyEDITH has heritage from the exposure time calculator built for the Altruistic Yield Optimizer (<a href='https://ui.adsabs.harvard.edu/abs/2014ApJ...795..122S/abstract'>C.C. Stark et al., 2014</a>), and has been validated against the AYO, exoSIMS, and EBS exposure time calculators." +
                    '<p><p align="justify">Uploaded spectra can be in either fixed width two-column ASCII (wave flux) or FITS format, where the spectrum is in the second HDU (HDU1, BINTABLE) with column 1 = "WAVELENGTH" and column 2 = "FLUX". Wavelength is assumed to be in Angstroms, Flux in FLAM (erg s−1 cm−2˚𝐴−1)</p>'
        )
        self.hrpanel1 = Div(text="<p>------------------- star ---------------------</p>")
        self.hrpanel2 = Div(text="<p>------------------ planet --------------------</p>")
        self.hrpanel3 = Div(text="<p>----------------------------------------------</p>")

        # currently unused, as we only have EAC1 working
        self.eac_buttons = RadioButtonGroup(labels=self.EACS, active=0)
        self.eac_buttons.on_change("active", self.eac_callback)

        self.newsnr  = Slider(title="Target SNR", value=10., start=0.1, end=100.0, step=0.1, ) 
        self.newsnr.on_change("value", self.snr_callback)

        self.newexp  = Slider(title="Target Exposure Time (hrs)", value=10, start=0.1, end=1000.0, step=0.1, )
        self.newexp.on_change("value", self.exp_callback)

        self.newdiameter  = Slider(title="Mirror Diameter", value=7., start=5, end=15, step=0.1, ) 
        self.newdiameter.on_change("value", self.diameter_callback)

        self.star = Select(title="Template Star Spectrum", value="G2V star", 
                    options=list(self.target_star.keys()), width=250) 
        self.star.on_change("value", self.star_callback)

        self.stellar_magnitude = Slider(title="V Magnitude of Star", value=12., start=-6, end=20, step=0.1, direction="rtl", sizing_mode="stretch_width") 
        self.stellar_magnitude.on_change("value", self.stellarmag_callback)

        self.stellar_radius = Slider(title="Radius of Star", value=1., start=.1, end=12., step=0.1, sizing_mode="stretch_width") 
        self.stellar_radius.on_change("value", self.stellar_radius_callback)

        self.distance  = Slider(title="Distance to System (pc)", value=10, start=1.4, end=100.0, step=0.1) 
        self.distance.on_change("value", self.distance_callback)

        self.planet = Select(title="Template Planet Spectrum", value="Earth", 
                options=list(self.target_planet.keys()), width=250)
        self.planet.on_change("value", self.planet_callback)

        self.semimajor = Slider(title="Semimajor Axis (AU)", value=0.1, start=0.01, end=10, step=0.01, ) 
        self.semimajor.on_change("value", self.semimajor_callback)

        self.delta_mag = Slider(title="delta Mag", value=15., start=10, end=30.0, step=0.1, ) 
        self.delta_mag.on_change("value", self.dmag_callback)

        self.upload = FileInput(accept=[".txt", ".csv", ".fit", ".fits", ".ascii", ".asdf"], title="Upload a Stellar Spectrum (.txt or FITS format, 10 MiB max)", directory=False, multiple=False) # 1. list allowed extensions
        self.upload.on_change("filename", self.process_spectrum)



    def tab_setup(self):
        self.observation_tab = TabPanel(child=self.exp_plot, title='Observation') # , width=400)
        self.info_panel = TabPanel(child=self.info_panel, title='Info') #, width=800)

        self.controls = column(children=[], sizing_mode='fixed', width=320, height=480) 
        self.starparam = row(children=[self.stellar_magnitude, self.stellar_radius], width_policy = "fit")

        self.exp_snr_toggle = RadioGroup(labels=["Solve For Exposure Time", "Solve For SNR"], active=0)
        def exp_snr_callback(active, old, new):
            if (new == 0):
                print(self.controls.children)
                self.controls.children = [self.intro, self.newdiameter, self.exp_snr_toggle, self.newsnr, self.hrpanel1, self.star, self.starparam, 
                                          self.distance, self.hrpanel2, self.planet, self.semimajor, self.hrpanel3, self.exptime_compute, self.upload, self.warning]
                outputs.tabs = [self.spec_panel, self.exp_panel, self.info_panel]
            elif new == 1:
                print(self.controls.children)
                self.controls.children = [self.intro, self.newdiameter, self.exp_snr_toggle, self.newexp, self.hrpanel1, self.star, self.starparam, 
                                          self.distance, self.hrpanel2, self.planet, self.semimajor, self.hrpanel3, self.snr_compute, self.upload, self.warning]
                outputs.tabs = [self.spec_panel, self.snr_panel, self.info_panel]                   
            #controls.change.emit()
            #outputs.change.emit()

        self.exp_snr_toggle.on_change("active", exp_snr_callback)

        # this is the initial for-exptime selection
        self.controls.children=[self.intro, self.newdiameter, self.exp_snr_toggle, self.newsnr, self.hrpanel1, self.star, self.starparam, self.distance, 
                                self.hrpanel2, self.planet, self.semimajor, self.hrpanel3, self.exptime_compute, self.upload, self.warning]

        outputs = Tabs(tabs=[self.spec_panel, self.exp_panel, self.info_panel], sizing_mode="inherit")
        plots = column(children=[outputs], sizing_mode='fixed', width=640, height=480)
        l = layout([[self.controls, plots]],sizing_mode='fixed', width=960, height=480)

        curdoc().theme = 'dark_minimal'
        curdoc().add_root(l) 
        curdoc().add_root(self.obsdata)

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

    def eac_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_eac": [self.EACS[new]], "observatory": [True]})

    def snr_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_snr": [new], "observation": [True]})
     
    def exp_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_exp": [new], "observation": [True]})

    def diameter_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_telescope_diameter": [new], "observatory": [True]})

    def star_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_star": [new], "scene": [True]})

    def stellarmag_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_stellar_magnitude": [new], "scene": [True]})

    def stellar_radius_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_stellar_radius": [new], "scene": [True]})

    def distance_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_distance": [new], "scene": [True]})

    def planet_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_planet": [new], "scene": [True]})

    def semimajor_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_semimajor": [new], "scene": [True]})

    def dmag_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_dMag": [new], "scene": [True]})

    def recalculate_exptime(self, newvalues):
        """
        The trick here is that Bokeh only synchronizes the calls at the end of a function, so
        if I want to change the button to "Please wait..." and THEN have it calculate, I have
        to make this call, and have this call fire off another callback with add_next_tick_callback.

        Parameters
        ----------
        newvalues : _type_
            _description_
        """
        self.exptime_compute.label = "Please Wait..."
        curdoc().add_next_tick_callback(partial(self.do_recalculate_exptime, newvalues))

    def recalculate_snr(self, newvalues):
        """
        The trick here is that Bokeh only synchronizes the calls at the end of a function, so
        if I want to change the button to "Please wait..." and THEN have it calculate, I have
        to make this call, and have this call fire off another callback with add_next_tick_callback.

        Parameters
        ----------
        newvalues : _type_
            _description_
        """

        self.snr_compute.label = "Please Wait..."
        curdoc().add_next_tick_callback(partial(self.do_recalculate_snr, newvalues))


    def process_spectrum(self, attr, old, new):
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
            #try:
                newstar = catalog.load_spec(f"../uploaded/{filename}", input_filename, filetype, magV=8, stellar_radius=1, planetary_radius=None, semimajor_axis=None, stargalaxy=True)
                if input_filename not in self.star.options:
                    self.star.options.append(input_filename)
                    self.target_star[input_filename] = newstar
                self.star.value = input_filename
                os.remove(f"../uploaded/{filename}") # don't clutter the upload directory
                self.star_callback([],[],input_filename)
            #except Exception as exc:
                #warning.text = str(exc)
        else:
            self.warning.text = "File too large"


coron_spec = CoronSpec()


coron_spec.load_initial()
