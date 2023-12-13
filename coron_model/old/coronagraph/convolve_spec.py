import numpy as np
import scipy as sp
from .degrade_spec import degrade_spec
from scipy import interp
from scipy import ndimage

__all__ = ['convolve_spec']

def convolve_spec(Ahr, lamhr, filters, forceTopHat=False):

    # if wl grid is backwards reverse it
    if lamhr[1] > lamhr[0]:
        pass
    else:
        lamhr=lamhr[::-1]
        Ahr=Ahr[::-1]

    # Sort filters by wavelength
    tdict = sorted(filters.__dict__.iteritems(), key=lambda x: x[1].bandcenter)
    F = []
    for x in tdict:
        if (x[1].wl == None) or (x[1].response == None) or forceTopHat:
            # Use FWHM with tophat convolution
            Fi = tophat_instrument(Ahr, lamhr, x[1].bandcenter, FWHM=x[1].FWHM)
        else:
            Fi = convolve_filter_response(lamhr, Ahr, x[1].wl, x[1].response, degrade=True)
        F.append(Fi)

    return np.array(F)

def convolve_filter_response(wlh, fh, wlf, response, degrade=False):

    # if wl grid is backwards reverse it
    if wlh[1] > wlh[0]:
        pass
    else:
        wlh=wlh[::-1]
        fh=fh[::-1]

    Nfilt = len(wlf)
    wlmin = np.min(wlf)
    wlmax = np.max(wlf)
    dlo = wlf[1] - wlf[0]
    dhi = wlf[-1] - wlf[-2]
    # Select only this wavelength region
    wlhr = wlh[(wlh > (wlmin-dlo*10)) & (wlh < (wlmax+dhi*10))]
    fhr = fh[(wlh > (wlmin-dlo*10)) & (wlh < (wlmax+dhi*10))]
    Nspec = len(wlhr)

    if degrade:
        try:
            F = degrade_spec(fhr, wlhr, wlf)
            R = response
        except ValueError:
            print 'Error in degrade_spec, switching to numpy interpolation.'
            if Nspec > Nfilt:
                # Interpolate filter response to hi-res spectral grid
                R = np.interp(wlhr, wlf, response)
                F = fhr
            else:
                # Interpolate spectrum to hi-res filter grid
                F = np.interp(wlf, wlhr, fhr)
                R = response
    elif Nspec > Nfilt:
        # Interpolate filter response to hi-res spectral grid
        R = np.interp(wlhr, wlf, response)
        F = fhr
    else:
        # Interpolate spectrum to hi-res filter grid
        F = np.interp(wlf, wlhr, fhr)
        R = response


    # Convolve with normalized filter response function
    F = F * (R / np.sum(R))

    # Integrate by summing
    Fc = np.sum(F)

    if Fc < 0.0:
        print Fc, wlf, degrade

    return Fc

def tophat_instrument(Fp, wl_hr, wlgrid, FWHM=0.035):

    Fratio11=tophatfold(wl_hr, Fp, FWHM)
    Fratio=interp(wlgrid,wl_hr,Fratio11)

    return Fratio

def tophatfold(lam, flux, FWHM=0.035):
    lammin=min(lam)
    lammax=max(lam)
    dlambda=FWHM/17.
    interlam=np.arange(lammin,lammax,dlambda)
    interflux=interp(interlam,lam,flux)

    #convovle flux array with gaussian--use smooth
    fold=sp.ndimage.filters.uniform_filter(interflux,size=17)

    #interpolate back to original grid
    fluxfold=interp(lam,interlam,fold)

    return fluxfold
