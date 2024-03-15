#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 15:04:24 2017

@author: gkanarek
"""
from __future__ import (print_function, division, absolute_import, with_statement,
                        nested_scopes, generators)

import astropy.units as u
#from astropy.modeling.functional_models import Box1D
#import numpy as np
import pysynphot as pys

#Define a new unit for spectral flux density
flambda = u.def_unit(["flambda","flam"], (u.photon / u.s / u.cm**2 / u.nm))

def renorm_sed(sed, new_magnitude, bandpass='johnson,v'):
    """
    Utility to renormalize an SED to a new manitude.
    """
    band = pys.ObsBandpass(bandpass)
    band.convert(sed.waveunits)
        
    new_sed = sed.renorm((new_magnitude + 2.5*u.mag('AB')).value, 'abmag', band)
    new_sed.convert('nm')
    new_sed.convert('abmag')
    
    return new_sed

def mag_from_sed(sed, camera):
    """
    Given a spectral energy distribution (SED) and a camera, convolve the
    SED with the camera's bandpass(es) to generate a set of AB magnitudes.
    
    NOTE AS OF 2017-11-20: not going to convolve or anything atm, just spit
    out pysynphot samples at the band pivotwaves. Convolution & integration
    over the camera bandpasses will be implemented in a future version.
    
    Parameters:
        sed    - pysynphot spectrum
        camera - Camera model
    """
    
    #Acquire camera bandpasses, making use of astropy.modeling model sets.
    pivots = camera.recover('pivotwave')
    
    sed.convert('ABMag')
    sed.convert(pivots.unit.name)
    output_mag = sed.sample(pivots.value)
    
    return output_mag * u.ABmag

    """## OLD VERSION:
    pivots = pivots.to(wavelength.unit)
    passes = passes.to(wavelength.unit)
    n = camera.n_bands
    bands = Box1D(amplitudes=[1.]*n, x_0=pivots.value, width=passes.value, n_models=n)
    
    #Calculate bandpass coefficients at the spectrum's wavelengths
    band_coeff = bands(np.tile(wavelength.value[None,...], (n,1)))
    
    #Convolve with flux and integrate
    convol = band_coeff * np.tile(sedflux.value[None,...], (n,1))
    integ = convol.sum(axis=0)
    
    #Ensure the proper unit conversions.
    integ = integ * sedflux.unit * wavelength.unit
    output = [i.to(u.ABmag, u.spectral_density(p)) for i,p in zip(integ, pivots)]
    output = np.array([o.value for o in output]) * u.ABmag
    
    return output"""
    