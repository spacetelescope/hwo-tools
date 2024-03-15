#!/usr/bin/env python
"""
Created on Tue Apr 25 11:44:58 2017

@author: gkanarek
"""

from __future__ import (print_function, division, absolute_import, with_statement,
                        nested_scopes, generators)

import astropy.units as u
import pysynphot as pys
import numpy as np

# JT omitted 050323 
#u.magnitude_zero_points.enable()

class JsonUnit(object):
    """
    A quick and dirty solution for making units storable via JSON. This is
    used in both the persistence subpackage (when implementing the JSON protocol),
    and in the actual tools, so that we can use astropy units along with
    Bokeh server.
    """
    
    def __init__(self, quant=0, unit=""):
        self._array = False
        if isinstance(quant, u.Quantity):
            if unit:
                quant = quant.to(u.Unit(unit))
            self._value, self._unit = self._grab_from_quantity(quant)
        else:
            self._unit = self._grab_unit(unit)
            self._value = self._grab_value(quant)
    
    def __repr__(self):
        if self._array:
            val, unit = self.value, self.unit
            return "{} {}".format(np.array_str(val, precision=3), unit)
        return str(self.use)

    def _grab_from_quantity(self, quantity):
        return self._grab_value(quantity.value), self._grab_unit(quantity.unit)
    
    def _grab_value(self, value):
        if isinstance(value, np.ndarray):
            self._array = True
            return value.tolist()
        else:
            self._array = False
            return value
    
    def _grab_unit(self, unit):
        if isinstance(unit, (u.UnitBase, u.FunctionUnitBase)):
            return unit.to_string()
        return unit
    
    @property
    def unit(self):
        return u.Unit(self._unit)
    
    @unit.setter
    def unit(self, new_unit):
        self._unit = self._grab_unit(new_unit)
    
    @property
    def value(self):
        if self._array:
            return np.array(self._value)
        return self._value
    
    @value.setter
    def value(self, new_val):
        self._value = self._grab_value(new_val)    
    
    @property
    def use(self):
        return self.value * self.unit
    
    @use.setter
    def use(self, new_quant):
        if not isinstance(new_quant, u.Quantity):
            raise TypeError("JsonUnit.use expects a Quantity object.")
        self._value, self._unit = self._grab_from_quantity(new_quant)
        
    def convert(self, new_unit):
        """
        Setting instance.unit changes the unit directly, and doesn't alter the
        value(s).
        """
    
        quant = self.use.to(new_unit)
        self._value, self._unit = self._grab_from_quantity(quant)
    
    def encode_json(self):
        return ("JsonUnit", {"unit": self._unit, "value": self._value})
    
    @classmethod
    def decode_json(cls, serialized):
        if isinstance(serialized, np.ndarray) or "JsonUnit" not in serialized:
            raise ValueError("Serialized element is not a JsonUnit")
        junit = cls()
        judict = serialized[1]
        junit._unit = judict["unit"]
        if isinstance(judict["value"],list):
            junit._array = True
        junit._value = judict["value"]
        return junit

class JsonSpectrum(object):
    """
    A version of JsonUnit to handle pysynphot spectra.
    """
    def __init__(self, spectrum=None):
        if spectrum is not None:
            self._wave = spectrum.wave
            self._wunit = spectrum.waveunits.name
            self._flux = spectrum.flux
            self._funit = spectrum.fluxunits.name
            
    
    def encode_json(self):
        return ['JsonSpectrum', {'wave': self._wave.tolist(),
                                 'wunit': self._wunit,
                                 'flux': self._flux.tolist(),
                                 'funit': self._funit}]
    
    @classmethod
    def decode_json(cls, serialized):
        if not isinstance(serialized, list) or "JsonSpectrum" not in serialized:
            raise ValueError("Serialized element is not a JsonSpectrum")
        jspec = cls()
        jsdict = serialized[1]
        jspec._wave = np.array(jsdict['wave'])
        jspec._wunit = jsdict['wunit']
        jspec._flux = np.array(jsdict['flux'])
        jspec._funit = jsdict['funit']
        return jspec
    
    @property
    def use(self):
        spec = pys.ArraySpectrum(wave=self._wave, flux=self._flux,
                                 waveunits=self._wunit, 
                                 fluxunits=self._funit)
        return spec
    
    @use.setter
    def use(self, new_spec):
        if not isinstance(new_spec, pys.spectrum.SourceSpectrum):
            raise TypeError("JsonSpectrum.use expects a pysynphot.spectrum.SourceSpectrum object.")
        self._wave = new_spec.wave
        self._wunit = new_spec.waveunits.name
        self._flux = new_spec.flux
        self._funit = new_spec.fluxunits.name
    
    def __repr__(self):
        wave, flux = self.wave, self.flux
        ww = "{} {}".format(np.array_str(wave.value, precision=3), wave.unit)
        ff = "{} {}".format(np.array_str(flux.value, precision=3), flux.unit)
        return "<Spectrum (wave {}; flux {})>".format(ww, ff)
    
    @property
    def wave(self):
        return self._wave * u.Unit(self._wunit)
    
    @wave.setter
    def wave(self, new_wave):
        if not isinstance(new_wave, list):
            new_wave = pre_encode(new_wave)
        if new_wave[0] == 'JsonUnit':
            self._wave = np.array(new_wave[1]['value'])
            if new_wave[1]['unit']:
                self._wunit = new_wave[1]['unit']
        else:
            self._wave = np.array(new_wave)
    
    @property
    def flux(self):
        if self._funit == 'abmag':
            return self._flux * u.ABmag
        return self._flux * u.Unit(self._funit)
    
    @flux.setter
    def flux(self, new_flux):
        if not isinstance(new_flux, list):
            new_flux = pre_encode(new_flux)
        if new_flux[0] == 'JsonUnit':
            self._flux = np.array(new_flux[1]['value'])
            if new_flux[1]['unit']:
                self._funit = new_flux[1]['unit']
        else:
            self._flux = np.array(new_flux)


def recover_quantities(*args):
    """
    Utility function to convert a number of JsonUnit instances into their
    Quantity counterparts, so that calculations may be done.
    """
    return [jsu.use for jsu in args]

def pre_encode(quant):
    """
    We only REALLY want to store the JSON serialization, a lot of the time.
    
    Same for JsonSpectra.
    """
    
    if isinstance(quant, pys.spectrum.SourceSpectrum):
        return JsonSpectrum(quant).encode_json()
    
    if not isinstance(quant, u.Quantity):
        return quant
    
    return JsonUnit(quant).encode_json()

def pre_decode(serialized):
    """
    Often, we don't actually need to store a JsonUnit representation of
    a quantity, we just want the Quantity version.
    
    Same for JsonSpectra.
    """
    try:
        spec = JsonSpectrum.decode_json(serialized)
        return spec.use
    except (AttributeError, ValueError, TypeError):
        pass #not a JsonSpectrum serialization
    
    try:
        quant = JsonUnit.decode_json(serialized)
    except (AttributeError, ValueError, TypeError):
        return serialized #not a JsonUnit serialization
    
    return quant.use

def str_jsunit(jsunit):
    try:
        spec = JsonSpectrum.decode_json(jsunit)
        return str(spec)
    except (AttributeError, ValueError, TypeError):
        pass #not a JsonSpectrum serialization
    
    try:
        quant = JsonUnit.decode_json(jsunit)
    except (AttributeError, ValueError, TypeError):
        return str(jsunit) #not a JsonUnit serialization
    
    return str(quant)
        
