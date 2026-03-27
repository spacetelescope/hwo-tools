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

from . import catalog

class pyEDITHETC():
    # classmethods
    target_planet, target_star = catalog.load_catalog()

    parameters = {}
    scene = pE.AstrophysicalScene()
    observation = pE.Observation() # define the observation object
    observatory = None # this piece, alone, has to be created WITH some configured parameters. So that's done in load_initial()
    obsdata = ColumnDataSource(data=dict(wavelength=[], exptime=[], FpFs=[], obs=[], noise_hi=[], noise_lo=[], snr=[]))
    inputs = ColumnDataSource(data=dict())

    def __init__(self):

        self.widget_setup()
        self.tab_setup()
    
    def widget_setup(self):

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

        self.stellar_diameter = Slider(title="Diameter of Star", value=1., start=.1, end=12., step=0.1, sizing_mode="stretch_width") 
        self.stellar_diameter.on_change("value", self.stellar_diameter_callback)

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
        self.starparam = row(children=[self.stellar_magnitude, self.stellar_diameter], width_policy = "fit")

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
        raise ValueError("Must be run from a customized class")

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
        self.inputs.data.update({"new_diameter": [new], "observatory": [True]})

    def star_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_star": [new], "scene": [True]})

    def stellarmag_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_stellar_magnitude": [new], "scene": [True]})

    def stellar_diameter_callback(self, attr, old, new):
        print(attr, old, new)
        self.inputs.data.update({"new_stellar_diameter": [new], "scene": [True]})

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

    def update_scene(self):
        #parameters.update(updates)

        self.scene.load_configuration(self.parameters)
        self.scene.calculate_zodi_exozodi(self.parameters)
        self.scene.validate_configuration()

    def update_observation(self):

        self.observation.load_configuration(self.parameters) # load the specified configuration in the parameters dict 
        self.observation.set_output_arrays()
        self.observation.validate_configuration()
        

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

    def recompute_planet_flux(self):
        print(self.parameters)
        solid_angle = self.parameters["planetary_radius"]**2/(4 * (self.parameters["semimajor_axis"]*1.5e8)**2) #Momentarily put both in km. pi cancels out of top and bottom. 
        flux_planet = self.parameters["FstarV_10pc"] * solid_angle * self.reflect_planet(self.parameters["wavelength"] << u.micron)
        self.parameters["F0"] = flux_planet.value
        self.parameters["Fp/Fs"] = (flux_planet / self.parameters["FstarV_10pc"]).value
        self.semimajor.value = self.parameters["semimajor_axis"] # Make sure it matches what was used


        print("Star:", self.parameters["FstarV_10pc"])
        print("Planet:", self.reflect_planet(self.parameters["wavelength"] << u.micron), solid_angle)
        print(self.parameters["semimajor_axis"], self.parameters["distance"])

    def load_planet(self, sourceID):
        # this is an albedo; the amount of incident flux received at that distance
        self.reflect_planet = self.target_planet[sourceID]["spectrum"]
        self.parameters["planetary_radius"] = self.target_planet[sourceID]["planetary_radius"]
        self.parameters["semimajor_axis"] = self.target_planet[sourceID]["semimajor_axis"]
        self.semimajor.value = self.target_planet[sourceID]["semimajor_axis"]

        print("Planet", self.parameters)
        

    def load_star(self, sourceID):
        star = self.target_star[sourceID]["spectrum"]
        self.parameters["magV"] = self.target_star[sourceID]["magV"] # Johnson V magnitude, specifically.
        self.parameters["stellar_radius"] = self.target_star[sourceID]["stellar_radius"]
        #magnitude.value = target_star[sourceID]["magV"]
        self.parameters["current_star"] = star


    def recompute_star_flux(self):
        bp = stsyn.band("johnson,v")

        self.stellar_magnitude.value = self.parameters["magV"]
        self.stellar_diameter.value = self.parameters["stellar_radius"]
        # we do not move the star from 10 pc, we merely provide the magnitude (and flux) at 10 pc and pyEDITH does the rest
        new_star = self.parameters["current_star"].normalize(self.parameters["magV"] * u.ABmag, band=bp, force="taper")
        flux = new_star(self.parameters["wavelength"]<< u.micron)
        #magnitude.value = parameters["magV"] # make sure it matches what was used

        self.parameters["Fstar_10pc"] = syn.units.convert_flux(self.parameters["wavelength"], flux, u.photon / (u.s * u.cm**2 * u.nm)).value
        # get the flux at 10 pc in the V band
        self.parameters["FstarV_10pc"] = syn.units.convert_flux(self.parameters["wavelength"], syn.Observation(new_star, bp, force="taper").effstim(), u.photon / (u.s * u.cm**2 * u.nm)).value

        #print("Star Flux", flux)

    # Ordinarily, these would be separate, but at the moment all changes here would seem to affect observatory, observation, and scene
    # it is, particularly, unclear what 
    def update_calculation(self, newvalues):
        print("------------------------------------")
        print(newvalues.data)

        if "new_star" in newvalues.data:
            print("Changed star")
            self.load_star(newvalues.data["new_star"][0])
            self.recompute_star_flux()
            self.recompute_planet_flux()
            del newvalues.data["new_star"] # consume the new value
        if "new_stellar_magnitude" in newvalues.data:
            print("Changed stellar magnitude")
            self.parameters["magV"] = newvalues.data["new_stellar_magnitude"][0]
            self.recompute_star_flux()
            del newvalues.data["new_stellar_magnitude"] # consume the new value
        if "new_stellar_diameter" in newvalues.data:
            print("Changed stellar diameter")
            self.parameters["stellar_radius"] = newvalues.data["new_stellar_diameter"][0]
            self.recompute_star_flux()
            del newvalues.data["new_stellar_diameter"] # consume the new value
        if "new_planet" in newvalues.data:
            print("Changed planet")
            self.load_planet(newvalues.data["new_planet"][0])
            self.recompute_planet_flux() # trigger a recomputation of the planetary flux
            del newvalues.data["new_planet"] # consume the new value
        if "new_semimajor" in newvalues.data:
            print("Changed Semimajor Axis")
            self.parameters["semimajor_axis"] = newvalues.data["new_semimajor"][0]
            self.recompute_planet_flux()
            del newvalues.data["new_semimajor"] # consume the new value
        if "new_distance" in newvalues.data:
            print("Changed System distance")
            self.parameters["distance"] = newvalues.data["new_distance"][0]
            self.recompute_planet_flux()
            del newvalues.data["new_distance"] # consume the new value
        if "new_snr" in newvalues.data:
            print("Changed SNR")
            self.parameters["snr"] = newvalues.data["new_snr"][0] * np.ones_like(self.parameters["wavelength"])
            del newvalues.data["new_snr"] # consume the new value
        if "new_exp" in newvalues.data:
            print("Changed Exposure Time")
            self.parameters["exptime"] = (newvalues.data["new_exp"][0] * np.ones_like(self.parameters["wavelength"]) << u.hr).to_value(u.s)
            del newvalues.data["new_exp"] # consume the new value
        if "new_eac" in newvalues.data:
            self.parameters["observatory_preset"] = newvalues.data["new_eac"][0]
            print("Changed EAC")
            del newvalues.data["new_eac"] # consume the new value
        if "new_telescope_diameter" in newvalues.data:
            print("Changed Telescope Diameter")
            self.parameters["diameter"] = newvalues.data["new_telescope_diameter"][0]
            del newvalues.data["new_diameter"] # consume the new value
        else:
            self.parameters["observatory_preset"] = "EAC1" # tells ETC to use EAC1 yaml files throughputs

        if "observation" in newvalues.data and newvalues.data["observation"][0]:
            print("Rerun observation...")
            self.update_observation()
            newvalues.data["observation"][0] = False

        if "scene" in newvalues.data and newvalues.data["scene"][0]:
            print("Rerun scene")
            self.update_scene()
            if self.parameters["regrid_wavelength"] is True:
                self.scene.regrid_spectra(self.parameters, self.observation)
            newvalues.data["scene"][0] = False

        pE.ObservatoryBuilder.configure_observatory(
            self.observatory, self.parameters, self.observation, self.scene
        )
        self.observatory.validate_configuration()
        #print(observatory.telescope.__dict__)
        #print_observatory(observatory)

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
