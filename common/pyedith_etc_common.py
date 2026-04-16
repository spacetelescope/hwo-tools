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

    def __init__(self):
        pass
    
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
        self.stellar_radius.value = self.parameters["stellar_radius"]
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
        if "new_stellar_radius" in newvalues.data:
            print("Changed stellar radius")
            self.parameters["stellar_radius"] = newvalues.data["new_stellar_radius"][0]
            self.recompute_star_flux()
            del newvalues.data["new_stellar_radius"] # consume the new value
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
