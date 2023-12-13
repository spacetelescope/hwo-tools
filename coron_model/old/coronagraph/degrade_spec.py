# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 11:00:41 2015

@author: jlustigy
"""

import numpy as np
import scipy as sp
from scipy import interpolate
from scipy.stats import binned_statistic

__all__ = ['degrade_spec', 'downbin_spec']

def degrade_spec(specHR, lamHR, lamLR, dlam=None):
    """
    Degrade high-resolution spectrum to lower resolution

    Parameters
    ----------
    specHR : array-like
        Spectrum to be degraded
    lamHR : array-like
        High-res wavelength grid
    lamLR : array-like
        Low-res wavelength grid
    dlam : array-like, optional
        Low-res wavelength width grid

    Returns
    -------
    specLO : ndarray
        Low-res spectrum
    """

    # Store input variables (not 100% necessary)
    lamHI  = np.array(lamHR)
    spec   = np.array(specHR)
    lamLO  = np.array(lamLR)
    if dlam is not None:
        dlamLO = dlam


    # Reverse ordering if wl vector is decreasing with index
    if lamHR[0] > lamHR[1]:
        lamHI = np.array(lamHR[::-1])
        spec = np.array(specHR[::-1])
    if lamLR[0] > lamLR[1]:
        lamLO = np.array(lamLR[::-1])
        if dlam is not None:
            dlamLO = dlam[::-1]

    # Number of gridpoints in output
    Nspec = len(lamLO)
    specLO = np.zeros(Nspec)

    # Loop over all spectral elements
    for i in range(Nspec):
        #print 'i=',i
        # If dlam not provided, must determine bin widths
        if dlam is None:
            # Define short and long wavelength edges of the
            # low res grid, special cases at the edge of the grid
            if i > 0 and i < Nspec-1:
                #print "loop1"
                lamS = 0.5*(lamLO[i-1] + lamLO[i])
                lamL = 0.5*(lamLO[i+1] + lamLO[i])
            elif i == 0:
                #print "loop2"
                lamS = lamLO[i] - 0.5*(lamLO[i+1] - lamLO[i])
                lamL = 0.5*(lamLO[i+1] + lamLO[i])
            else:
                #print "loop3"
                lamS = 0.5*(lamLO[i-1] + lamLO[i])
                lamL = lamLO[i] + 0.5*(lamLO[i] - lamLO[i-1])
        else:
            lamS = lamLO[i] - 0.5*dlamLO[i]
            lamL = lamLO[i] + 0.5*dlamLO[i]

        # Elements of hi-res grid within the low-res element
        eps = 1e-10
        iss = (lamHI - lamS >= eps) & (lamHI - lamL <= eps)

        # If there aren't any high-res elements within low
        # res element, then error
        check1 = False if True in iss else True
        if check1:
            print("Error in DEGRADE_SPECTRUM: wavelength grids do not sync")

        # If 3 of less elements of spectrum within low-res gridpoint,
        # then do an interpolation, otherwise, integrate the high-res
        # spectrum over the low-res element, distributing the integrated
        # energey into the low-res element
        if len(lamHI[iss]) == 0:
            print "No HiRes elements in Low Res bin!"
            import sys; sys.exit()
        elif len(lamHI[iss]) == 1:
            specs = spec[iss]
        elif len(lamHI[iss]) <= 3:
            interpfunc = interpolate.interp1d(lamHI[iss], spec[iss], kind='linear',bounds_error=False)
            specs = interpfunc(lamLO[i])
        else:
            interpfunc = interpolate.interp1d(lamHI[iss], spec[iss], kind='linear',bounds_error=False, fill_value=0.0)
            speci = np.hstack([interpfunc(lamS), spec[iss], interpfunc(lamL)])
            #print speci
            lami = np.hstack([lamS,lamHI[iss],lamL])
            #print lami
            specs = np.trapz(speci,x=lami) / (lamL - lamS)
            #print specs

        # Insert result into output array
        specLO[i] = specs

    return specLO

def downbin_spec(specHR, lamHR, lamLR, dlam=None):
    """
    Re-bin spectum to lower resolution using scipy.binned_statistic

    Parameters
    ----------
    specHR : array-like
        Spectrum to be degraded
    lamHR : array-like
        High-res wavelength grid
    lamLR : array-like
        Low-res wavelength grid
    dlam : array-like, optional
        Low-res wavelength width grid

    Returns
    -------
    specLR : ndarray
        Low-res spectrum
    """

    if dlam is None:
        ValueError("Please supply dlam in downbin_spec()")

    # Reverse ordering if wl vector is decreasing with index
    if lamHR[0] > lamHR[1]:
        lamHI = np.array(lamHR[::-1])
        spec = np.array(specHR[::-1])
    if lamLR[0] > lamLR[1]:
        lamLO = np.array(lamLR[::-1])
        dlamLO = np.array(dlam[::-1])

    # Calculate bin edges
    LRedges = np.hstack([lamLR - 0.5*dlam, lamLR[-1]+0.5*dlam[-1]])

    # Call scipy.stats.binned_statistic()
    specLR = binned_statistic(lamHR, specHR, statistic="mean", bins=LRedges)[0]

    return specLR
