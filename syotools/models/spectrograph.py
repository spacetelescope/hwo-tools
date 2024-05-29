#!/usr/bin/env python
"""
Created on Sat Oct 15 16:56:40 2016

@author: gkanarek, tumlinson
"""

from __future__ import (print_function, division, absolute_import, with_statement,
                        nested_scopes, generators)

import numpy as np
import astropy.units as u
from astropy.table import QTable

from syotools.models.base import PersistentModel
from syotools.models.exposure import SpectrographicExposure
from syotools.defaults import default_spectrograph
from syotools.defaults import default_spectropolarimeter
from syotools.utils import pre_encode

class Spectrograph(PersistentModel):
    """
    The basic spectrograph class, which provides parameter storage for 
    optimization.
    
    Attributes: #adapted from the original in Telescope.py
        telescope    - the Telescope object associated with this spectrograph
        exposures    - the list of Exposures taken with this spectrograph
    
        name         - name of the spectrograph (string)
        
        modes        - supported observing modes (list)
        descriptions - description of supported observing modes (dict)
        mode         - current observing mode (string)
        bef          - background emission function in erg/cm2/pixel/s ement (float array)
        R            - spectral resolution (float)
        wrange        - effective wavelength range (2-element float array)
        wave         - wavelength in Angstroms (float array)
        aeff         - effective area at given wavelengths in cm^2 (float array)
        
        _lumos_default_file - file path to the fits file containing LUMOS values
        
        _default_model - used by PersistentModel
    """
    
    _default_model = default_spectrograph
    
    telescope = None
    exposures = []
    
    _lumos_default_file = ''
    
    name = ''
    modes = []
    descriptions = {}
    # JT fixed BEF units 05292024 
    #bef = pre_encode(np.zeros(0, dtype=float) * (u.erg / u.cm**2 / u.s / u.AA))
    bef = pre_encode(np.zeros(0, dtype=float) * (u.erg / u.cm**2 / u.s / u.pix))
    R = pre_encode(0. * u.dimensionless_unscaled) 
    wave = pre_encode(np.zeros(0, dtype=float) * u.AA)
    aeff = pre_encode(np.zeros(0, dtype=float) * u.cm**2)
    wrange = pre_encode(np.zeros(2, dtype=float) * u.AA)
    _mode = ''
    
    #Property wrapper for mode, so that we can use a custom setter to propagate
    #mode updates to all the rest of the parameters
    
    @property
    def mode(self):
        return self._mode
    
    @mode.setter
    def mode(self, new_mode):
        """
        Mode is used to set all the other parameters
        """
        
        nmode = new_mode.upper()
        if self._mode == nmode or nmode not in self.modes:
            return
        self._mode = nmode
        table = QTable.read(self._lumos_default_file, nmode)
                
        self.R = pre_encode(table.meta['R'] * u.pix)
        self.wave = pre_encode(table['Wavelength'])
        #JT 05292024 
        #self.bef = pre_encode(table['BEF'] / self.recover('delta_lambda'))
        self.bef = pre_encode(table['BEF']) 
        self.aeff = pre_encode(table['A_Eff'])
        wrange = np.array([table.meta['WAVE_LO'], table.meta['WAVE_HI']]) * u.AA
        self.wrange = pre_encode(wrange)

    @property
    def delta_lambda(self):
        wave, R = self.recover('wave', 'R')
        return pre_encode(wave / R)
    
    def create_exposure(self):
        new_exposure = SpectrographicExposure()
        self.add_exposure(new_exposure)
        return new_exposure
    
    def add_exposure(self, exposure):
        self.exposures.append(exposure)
        exposure.spectrograph = self
        exposure.telescope = self.telescope
        exposure.calculate()

class Spectropolarimeter(Spectrograph):
    """
    The basic spectropolarimeter class for POLLUX, which provides parameter storage for 
    optimization.
    
    Attributes: #adapted from the original in Telescope.py
        telescope    - the Telescope object associated with this spectrograph
        exposures    - the list of Exposures taken with this spectrograph
    
        name         - name of the spectrograph (string)
        
        modes        - supported observing modes (list)
        descriptions - description of supported observing modes (dict)
        mode         - current observing mode (string)
        bef          - background emission function in erg/cm2/s/pixel (float array)
        R            - spectral resolution (float)
        wrange        - effective wavelength range (2-element float array)
        wave         - wavelength in Angstroms (float array)
        aeff         - effective area at given wavelengths in cm^2 (float array)
        
        _lumos_default_file - file path to the fits file containing LUMOS values
        
        _default_model - used by PersistentModel
    """
    
    _default_model = default_spectropolarimeter

