"""
Created on Mon Jul  6 11:00:41 2015

@author: jlustigy
"""

import numpy as np
import scipy as sp
from scipy import interpolate

def degrade_spec(specHR, lamHR, lamLR, dlam=[]):
    
    # Store input variables (not 100% necessary) 
    lamHI  = np.array(lamHR)
    spec   = np.array(specHR)
    lamLO  = np.array(lamLR)
    if len(dlam) != 0:
        dlamLO = dlam 
    
    
    # Reverse ordering if wl vector is decreasing with index
    if lamHR[0] > lamHR[1]:
        lamHI = np.array(lamHR[::-1])
        spec = np.array(specHR[::-1])
    if lamLR[0] > lamLR[1]:
        lamLO = np.array(lamLR[::-1])
        if len(dlam) != 0:
            dlamLO = dlam[::-1] 
    
    # Number of gridpoints in output
    Nspec = len(lamLO)
    specLO = np.zeros(Nspec)
    
    # Loop over all spectral elements
    for i in range(Nspec):
        #print 'i=',i
        # If dlam not provided, must determine bin widths
        if len(dlam) == 0:
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
        if len(lamHI[iss]) <= 3:
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
























