#!/usr/bin/env python
"""
Created on Fri Oct 14 21:31:18 2016

@author: gkanarek, tumlinson
"""

from __future__ import (print_function, division, absolute_import, with_statement,
                        nested_scopes, generators)

import numpy as np
import astropy.constants as const
import astropy.units as u

from syotools.models.base import PersistentModel
from syotools.models.exposure import PhotometricExposure
from syotools.defaults import default_camera
from syotools.utils import pre_encode, pre_decode
from syotools.spectra.utils import mag_from_sed

def nice_print(arr):
    """
    Utility to make the verbose output more readable.
    """
    
    arr = pre_decode(arr) #in case it's a JsonUnit serialization

    if isinstance(arr, u.Quantity):
        l = ['{:.2f}'.format(i) for i in arr.value]
    else:
        l = ['{:.2f}'.format(i) for i in arr]
    return ', '.join(l)

class Camera(PersistentModel): 
    """
    The basic camera class, which provides parameter storage for 
    optimization.
    
    Attributes: #adapted from the original in Telescope.py
        telescope    - the Telescope object associated with the camera
        exposures    - the list of Exposures taken with this camera
    
        name         - name of the camera (string)
        n_bands      - number of wavelength bands (int)
        n_channels   - number of channels (int)
        pivotwave    - central wavelengths for bands, in nanometers (float array)
        bandnames    - names of bands (string list)
        channels     - grouping of bands into channels [UV, Optical, IR], 
                       and indicating the reference band for pixel size (list of tuples)
        fiducials    - fiducial wavelength of the band, for reference (float array)
        total_qe     - total quantum efficiency in each band (float array)
        ap_corr      - magnitude correction factor for aperture size (float array)
        bandpass_r   - resolution in each bandpass (float array)
        dark_current - dark current values in each band (float array)
        detector_rn  - read noise for the detector in each band (float array)
        sky_sigma    - sky background emission (float array)
        
        _default_model - used by PersistentModel
        
    The following are attributes I haven't included, and the justification:
        R_effective - this doesn't seem to be used anywhere
    """
    
    _default_model = default_camera
    
    telescope = None
    exposures = []
    
    name = ''
    pivotwave = pre_encode(np.zeros(1, dtype=float) * u.nm)
    bandnames = ['']
    channels = [([],0)]
    fiducials = pre_encode(np.zeros(1, dtype=float) * u.nm)
    total_qe = pre_encode(np.zeros(1, dtype=float) * u.dimensionless_unscaled)
    ap_corr = pre_encode(np.zeros(1, dtype=float) * u.dimensionless_unscaled)
    bandpass_r = pre_encode(np.zeros(1, dtype=float) * u.dimensionless_unscaled)
    dark_current = pre_encode(np.zeros(1, dtype=float) * (u.electron / u.s / u.pixel))
    detector_rn = pre_encode(np.zeros(1, dtype=float) * (u.electron / u.pixel)**0.5)
    sky_sigma = pre_encode(np.zeros(1, dtype=float) * u.dimensionless_unscaled)
    
    @property
    def pixel_size(self):
        """
        Compute the pixel size as a function of pivot wavelength.
        
        Use the reference band for each channel as the fiducial: U-band for UV 
        and optical, J-band for IR.
        """
        
        pixsize = np.zeros(self.n_bands, dtype=float)
        
        #Convert from JsonUnit to Quantity for calculation purposes.
        fiducials, aperture = self.recover('fiducials', 'telescope.aperture')
        
        for ref, bands in enumerate(self.channels):
            pxs = (0.5 * fiducials[ref] * u.rad / aperture).to(u.arcsec).value
            pixsize[bands] = pxs
        
        #serialize with JsonUnit for transportation purposes.
        return pre_encode(pixsize * u.arcsec / u.pix)
    
    @property
    def n_bands(self):
        return len(self.bandnames)
    
    @property
    def n_channels(self):
        return len(self.channels)
    
    @property
    def derived_bandpass(self):
        """
        Calculate the bandpasses.
        """
        
        #Convert to Quantity for calculations.
        pivotwave, bandpass_r = self.recover('pivotwave','bandpass_r')
        
        #serialize with JsonUnit for transportation.
        return pre_encode(pivotwave / bandpass_r)
    
    @property
    def ab_zeropoint(self):
        """
        AB-magnitude zero points as per Marc Postman's equation.
        """
        pivotwave = self.recover('pivotwave').to(u.nm)
        abzp = 5509900. * (u.photon / u.s / u.cm**2) / pivotwave
        
        return pre_encode(abzp)
        
    
    @property
    def fwhm_psf(self):
        """
        Calculate the FWHM of the camera's PSF.
        """
        #Convert to Quantity for calculations.
        pivotwave, aperture = self.recover('pivotwave', 'telescope.aperture')
        diff_limit, diff_fwhm = self.recover('telescope.diff_limit_wavelength',
                                             'telescope.diff_limit_fwhm')
        
        #fwhm = (1.22 * u.rad * pivotwave / aperture).to(u.arcsec)
        fwhm = (1.03 * u.rad * pivotwave / aperture).to(u.arcsec)
        
        #only use these values where the wavelength is greater than the diffraction limit
        fwhm = np.where(pivotwave > diff_limit, fwhm, diff_fwhm) * u.arcsec
        
        #serialize with JsonUnit for transportation.
        return pre_encode(fwhm)
    
    def _print_initcon(self, verbose):
        if verbose: #These are our initial conditions
            print('Telescope diffraction limit: {}'.format(pre_decode(self.telescope.diff_limit_wavelength)))
            print('Telescope aperture: {}'.format(pre_decode(self.telescope.aperture)))
            print('Telescope temperature: {}'.format(pre_decode(self.telescope.temperature)))
            print('Pivot waves: {}'.format(nice_print(self.pivotwave)))
            print('Pixel sizes: {}'.format(nice_print(self.pixel_size)))
            print('AB mag zero points: {}'.format(nice_print(self.ab_zeropoint)))
            print('Quantum efficiency: {}'.format(nice_print(self.total_qe)))
            print('Aperture correction: {}'.format(nice_print(self.ap_corr)))
            print('Bandpass resolution: {}'.format(nice_print(self.bandpass_r)))
            print('Derived_bandpass: {}'.format(nice_print(self.derived_bandpass)))
            print('Detector read noise: {}'.format(nice_print(self.detector_rn)))
            print('Dark rate: {}'.format(nice_print(self.dark_current)))
            
    def _fsky(self, verbose=True):
        """
        Calculate the sky flux as per Eq 6 in the SNR equation paper.
        """
        
        (f0, D, dlam, Phi, fwhm, Sigma) = self.recover('ab_zeropoint', 
                'telescope.effective_aperture', 'derived_bandpass', 
                'pixel_size', 'fwhm_psf', 'sky_sigma')
        
        D = D.to(u.cm)
        m = 10.**(-0.4 * Sigma) / u.arcsec**2
        Npix = self._sn_box(False)
        
        if verbose:
            print('Sky brightness: {}'.format(nice_print(Sigma)))
        
        fsky = f0 * np.pi / 4. * D**2 * dlam * m * (Phi**2 * Npix) * u.pix
        
        return fsky
    
    def _sn_box(self, verbose):
        """
        Calculate the number of pixels in the SNR computation box.
        """

        (Phi, fwhm_psf) = self.recover('pixel_size', 'fwhm_psf')
        sn_box = np.round(3. * fwhm_psf / Phi)
        
        if verbose:
            print('PSF width: {}'.format(nice_print(fwhm_psf)))
            print('SN box width: {}'.format(nice_print(sn_box)))
        
        return sn_box**2 / u.arcsec / u.arcsec / u.pix #don't want pix**2 units
    
    def c_thermal(self, verbose=True):
        """
        Calculate the thermal emission counts for the telescope.
        """
        
        #Convert to Quantities for calculation.
        (bandpass, pivotwave, aperture, ota_emissivity, 
         total_qe, pixel_size) = self.recover('derived_bandpass', 'pivotwave', 
                'telescope.effective_aperture',  'telescope.ota_emissivity', 
                'total_qe', 'pixel_size')
        
        box = self._sn_box(verbose)
        
        bandwidth = bandpass.to(u.cm)
    
        h = const.h.to(u.erg * u.s) # Planck's constant erg s 
        c = const.c.to(u.cm / u.s) # speed of light [cm / s] 
    
        energy_per_photon = h * c / pivotwave.to(u.cm) / u.ph
    
        D = aperture.to(u.cm) # telescope diameter in cm 
        
        Omega = (pixel_size**2 * box * u.pix).to(u.sr)
        
        planck = pre_decode(self.planck)
        qepephot = total_qe * planck / energy_per_photon
        
        if verbose:
            print('Planck spectrum: {}'.format(nice_print(planck)))
            print('QE * Planck / E_phot: {}'.format(nice_print(qepephot)))
            print('E_phot: {}'.format(nice_print(energy_per_photon)))
            print('Omega: {}'.format(nice_print(Omega)))
    
        thermal = (ota_emissivity * planck / energy_per_photon * 
    			(np.pi / 4. * D**2) * total_qe * Omega * bandwidth )
        
        #serialize with JsonUnit for transportation
        return pre_encode(thermal) 
    
    @property
    def planck(self):
        """
        Planck spectrum for the various wave bands.
        """
        #Convert to Quantities for calculation
        pivotwave, temperature = self.recover('pivotwave', 'telescope.temperature')
        
        wave = pivotwave.to('cm')
        temp = temperature.to('K')
        h = const.h.to(u.erg * u.s) # Planck's constant erg s 
        c = const.c.to(u.cm / u.s) # speed of light [cm / s] 
        k = const.k_B.to(u.erg / u.K) # Boltzmann's constant [erg deg K^-1] 
        x = 2. * h * c**2 / wave**5 
        exponent = (h * c / (wave * k * temp)) 
    
        result = (x / (np.exp(exponent)-1.)).to(u.erg / u.s / u.cm**3) / u.sr
        
        #serialize with JsonUnit for transportation
        return pre_encode(result)
    
    def interpolate_at_bands(self, sed):
        """
        Interpolate an SED to obtain magnitudes for the camera's wavebands.
        """
        return mag_from_sed(sed, self)
    
    def create_exposure(self):
        new_exposure = PhotometricExposure()
        self.add_exposure(new_exposure)
        return new_exposure
        
    def add_exposure(self, exposure):
        self.exposures.append(exposure)
        exposure.camera = self
        exposure.telescope = self.telescope
        exposure.calculate()
            
