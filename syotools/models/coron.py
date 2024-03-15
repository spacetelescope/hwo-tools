#!/usr/bin/env python
"""
Created on Tue Oct 18 15:23:49 2016

@author: gkanarek
"""

from __future__ import (print_function, division, absolute_import, with_statement,
                        nested_scopes, generators)

import astropy.units as u
import numpy as np

import syotools.coronagraph as cg
from syotools.models.base import PersistentModel
from syotools.defaults import default_coronagraph
from syotools.utils import pre_encode

class Coronagraph(PersistentModel):
    """
    The basic coronagraph class, which provides parameter storage for 
    optimization.
    
    Unlike other models, most of the calculations for this will be handled by
    the coronagraph package. This is mostly for parameter storage and initial
    conditions prep.
    
    Attributes: #adapted from the original in coron_model.py
        engine       - coronograph module instance
        count_rates  - output dict from coronograph.count_rates
    
        int_time     - Integration time for Bokeh slider, in hours (float)
        phase_angle  - planetary phase angle at quadrature, in degrees (float)
        phase_func   - planetary phase function at quadrature (already included in SMART run) (float)
        r_planet     - planetary radius in Earth radii for Bokeh slider (float)
        semimajor    - semi-major orbital axis in AU for Bokeh slider (float)
        t_eff        - Stellar effective temperature in K for Sun-like star (float)
        r_star       - Stellar radius in solar radii (float)
        d_system     - distance to planetary system in pc for Bokeh slider (float)
        n_exoz       - number of exo-zodis for Bokeh slider (float)
        wave         - hi-resolution wavelength array for Earth spectrum in microns (float array)
        radiance     - hi-resolution radiance array for Earth spectrum in W/m**2/um/sr (float array)
        sol_flux     - hi-resolution solar flux array for Earth spectrim in W/m**2/um (float array)
        
        _default_model - used by PersistentModel
    """
    
    _default_model = default_coronagraph
    
    engine = cg
    count_rates = {}
    
    int_time = pre_encode(0. * u.hr)
    phase_angle = pre_encode(0. * u.deg)
    phase_func = pre_encode(1.)
    r_planet = pre_encode(1. * u.R_earth)
    semimajor = pre_encode(1. * u.au)
    t_eff = pre_encode(0. * u.K)
    r_star = pre_encode(1. * u.R_sun)
    d_system = pre_encode(0. * u.pc)
    n_exoz = 1.
    wave = pre_encode(np.zeros(0, dtype=float) * u.um)
    radiance = pre_encode(np.zeros(0, dtype=float) * (u.W / u.m**2 / u.um / u.sr))
    sol_flux = pre_encode(np.zeros(0, dtype=float) * (u.W / u.m**2 / u.um))
    
    
    def __init__(self, *arg, **kw):
        super().__init__(*arg, **kw)
        self._calc_count_rates()
    
    @property
    def albedo(self):
        """
        Planetary albedo spectrum.
        """
        return pre_encode(np.pi * u.sr * (np.pi * self.radiance / self.sol_flux).decompose())
    
    def _calc_count_rates(self):
        """
        Compute the coronagraphic model using the coronagraph package.
        """
        #Convert to Quantities for calculations.
        (albedo, wave, sol_flux, phase_angle, phase_func, 
         r_planet, t_eff, r_star, semimajor, d_system) = self.recover('albedo',
                'wave', 'sol_flux', 'phase_angle', 'phase_func', 'r_planet',
                't_eff', 'r_star', 'semimajor', 'd_system')
        
        al = albedo.value
        wv = wave.to(u.um).value
        sf = sol_flux.to(u.W / u.m**2 / u.um).value
        pa = phase_angle.to(u.deg).value
        pf = phase_func.value
        rp = r_planet.to(u.R_earth).value
        te = t_eff.to(u.K).value
        rs = r_star.to(u.R_sun).value
        sm = semimajor.to(u.au).value
        ds = d_system.to(u.pc).value
        ez = self.n_exoz
        cr = cg.count_rates(al, wv, sf, pa, pf, rp, te, rs, sm, ds, ez)
        self._count_rates = dict(*zip(['wavelength','wave_bin','albedo',
                                       'quant_eff','flux_ratio','planet_cr',
                                       'speckle_cr','zodi_cr','exoz_cr',
                                       'dark_cr','read_cr','thermal_cr',
                                       'dtsnr'], cr))
    
    @property
    def background_cr(self):
        """
        Background photon count rates.
        """
        return sum([self._count_rates[x] for x in ['zodi_cr', 'exoz_cr', 
                                                   'speckle_cr', 'dark_cr',
                                                   'read_cr', 'thermal_cr']])
    
    @property
    def dts(self):
        """
        Integration time in seconds
        """
        return self.int_time.to(u.s).value
    
    @property
    def snr(self):
        """
        Calculate the SNR based on count rates.
        """
        pcr = self._count_rates['planet_cr']
        bcr = self.background_cr
        dts = self.dts
        return pcr * dts / np.sqrt((pcr + 2*bcr) * dts)
    
    @property
    def sigma(self):
        """
        Calculate the 1-sigma errors.
        """
        return self._count_rates['flux_ratio'] / self.snr

    @property
    def spectrum(self):
        """
        Create a spectrum by adding random noise to the flux ratio.
        """
        c_ratio = self._count_rates['flux_ratio']
        return c_ratio + np.random.randn(len(c_ratio)) * self.sigma

    @property
    def planet(self):
        """
        Generate the planet data dictionary for Bokeh.
        """
        spec = self.spectrum
        sig = self.sigma
        return {'lam':self._count_rates['wavelength'],
                'cratio': self._count_rates['flux_ratio'] * 1.e9,
                'spec': spec * 1.e9,
                'downerr': (spec - sig) * 1.e9,
                'uperr': (spec + sig) * 1.e9}