#!/usr/bin/env python
"""
Created on Tue Oct 18 15:23:49 2016
Updated Sept 2024 JT 

@author: gkanarek
"""

from __future__ import (print_function, division, absolute_import, with_statement,
                        nested_scopes, generators)

import astropy.units as u
import numpy as np

import syotools.coronagraph as ty
from syotools.models.base import PersistentModel
from syotools.models.exposure import CoronagraphicExposure
from syotools.defaults import default_coronagraph
from syotools.utils import pre_encode

class Coronagraph(PersistentModel):

# significant changes needed here to align with Mennesson paper !! 
# JT 9/10/2024 

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
    
    telescope = None 
    exposures = []

    engine = ty
    count_rates = {}
    
    int_time = pre_encode(10. * u.hr)
    int_time = 10. * u.hr 
    raw_contrast = 1e-10 # dimensionless 
    eta_p = 0.2  # core throughput, taking the value for a type 1 coronagraph from Mennesson Figure 24 
    Tr = 0.3 # total througput absent masks and detector QE - used number for all coronagraphs per Table 2
    quantum_eff = 0.9 # detector quantum efficiency (0.9 per Table 2) 

    sigma_DeltaC = 5e-12 # the post-processing factor 

    # all this below was needed only for calling the Ty code 
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
    #end of stuff for Ty code 
    
    
    def __init__(self, *arg, **kw):
        super().__init__(*arg, **kw)
        #self._calc_count_rates_imaging() 
        #self._calc_count_rates_ty() # Gray's code - commented out bc of error and in favor of new imaging mode  
    
    @property
    def albedo(self):
        """
        Planetary albedo spectrum.
        """
        return pre_encode(np.pi * u.sr * (np.pi * self.radiance / self.sol_flux).decompose())
    
    def _calc_count_rates_imaging(self):
        """
        Compute the coronagraphic model for a corongraph image using the math from Mennesson et al. (2024) 
        For now we are ignoring the xy dependence of these calculations, making this a zero-dimensional calculation.  
        """
        epsilon = 1.0e-10 # mean astrophysical planet-to-star flux ratio over the bandpass 

        # Sun at 10 pc in V band (V_abs_Sun = 4.8) in units photons cm^-2 s^-1 micron^-1 - using https://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/pet/magtojy/ 
        sun_flux = 4.51e-7 * u.erg / u.cm**2 / u.s / u.micron
        ph_energy = 6.626196e-27 * u.erg * u.s * (2.9979e10 * u.cm / u.s) / (5500 * u.Angstrom).to(u.cm) / u.ph 
        sun_phot = sun_flux / ph_energy 
 
        phi_S = sun_phot / 12.**2 * 10.**2 # correct photon rate to 12 pc for a Solar twin - should be photons / s/ cm^2 / micron 

        collecting_area = np.pi * (self.telescope.aperture.to(u.cm) / 2.)**2 # square cm, 6 meters, filled aperture, this will have to come in from the telescope object in the end. 

        dlambda = 0.55 / 5 * u.micron # bandpass taken from 0.55 micron central wavelength and spectral resolution per Table 2 
        
        N_s = phi_S * collecting_area * self.quantum_eff * self.Tr * dlambda # total number of stellar photons per sec in the absence of any corongraph - Eq (8) 

        # the various count rates 
        r_pl = self.eta_p * epsilon * N_s # planet photon rate (ph / s) detected in a circular photometric aperture Eq (3) 
        r_sp = self.eta_p * self.raw_contrast * N_s      # photon rate from residual starlight in speckles Eq (4) 
        r_sz = 1. * u.ph / u.s / u.micron * dlambda * self.Tr # read off the upper right panel of Figure 23. 
        r_xz = 7. * u.ph / u.s / u.micron * dlambda * self.Tr # read off the upper right panel of Figure 23. 

        r_n = r_pl + r_sp + r_sz + r_xz # total noise rate, sum of terms just below Eq (2) 
        post_term = self.eta_p * self.sigma_DeltaC * N_s * self.int_time.to(u.s) / u.ph #<---- see Eq 13

        # the factor of 2 on r_n here in the denominator is from Eq 13 and accounts for differential imaging 
        self._signal_to_noise = (r_pl * self.int_time.to(u.s) / u.ph) / (2. * r_n * self.int_time.to(u.s)/u.ph  + post_term**2 )**0.5 

    def _calc_count_rates_ty(self):
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
        # this is where the Ty & Giada code is called 
        cr = ty.count_rates(al, wv, sf, pa, pf, rp, te, rs, sm, ds, ez)
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

    def create_exposure(self):
        new_exposure = CoronagraphicExposure()
        self.add_exposure(new_exposure)
        return new_exposure

    def add_exposure(self, exposure):
        self.exposures.append(exposure)
        exposure.coronagraph = self
        exposure.telescope = self.telescope
        exposure.calculate()
