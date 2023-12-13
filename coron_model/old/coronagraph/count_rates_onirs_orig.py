# Import dependent modules
import numpy as np
import sys
from .degrade_spec import degrade_spec, downbin_spec
from .convolve_spec import convolve_spec
from .noise_routines import Fobj, cobj, czodi,  \
    cdark, cread, ctherm, ccic, f_airy, ctherm_earth, construct_lam, \
    set_quantum_efficiency, set_read_noise, set_dark_current, set_lenslet, \
    set_throughput, set_atmos_throughput, get_thermal_ground_intensity, \
    exptime_element, lambertPhaseFunction
import pdb
import os

__all__ = ['count_rates_onirs']

def count_rates_onirs(Fohr, lamhr,
                Ro, d,
                mode   = "IFS",
                filter_wheel = None,
                Phi = -1, #chosen as a value a user would never set 
                lammin = 0.4,
                lammax = 2.5,
                Res    = 70.0,
                diam   = 10.0,
                Tput   = 0.15,
                C      = 1e-10,
                IWA    = 3.0,
                OWA    = 20.0,
                Tsys   = 150.0,
                Tdet   = 50.0,
                emis   = 0.09,
                De     = 1e-4,
                DNHpix = 3.0,
                Re     = 0.1,
                Dtmax  = 1.0,
                X      = 1.5,
                qe     = 0.9,
                MzV    = 23.0, #23
                MezV   = 22.0, #22
                Res_NIR = -1,  #chosen as a value a user would never set 
                Res_UV = -1,  #chosen as a value a user would never set 
                wantsnr=10.0, FIX_OWA = False, COMPUTE_LAM = False,
                SILENT = False, NIR = True, UV=True, THERMAL = True,
                GROUND = False):
    """
    Runs coronagraph model (Robinson et al., 2016) to calculate planet and noise
    photon count rates for specified telescope and system parameters.

    Parameters
    ----------
    Fohr : array
        High-res, wavelength-dependent object flux array [W/m**2/um]
    lamhr : array
        High-res wavelength grid  [um]
    Ro : float
        Object radius [R_sun]
    d : float
        Distance to observed star-planet system [pc]
    mode : str, optional
        Telescope observing mode: "IFS" or "Imaging"
    filter_wheel : Wheel, optional
        Wheel object containing imaging filters
    lammin : float, optional
        Minimum wavelength [um]
    lammax : float, optional
        Maximum wavelength [um]
    Res : float, optional
        Instrument spectral resolution (lambda / delta_lambda)
    diam : float, optional
        Telescope diameter [m]
    Tput : float, optional
        Telescope and instrument throughput
    C : float, optional
        Coronagraph design contrast
    IWA : float, optional
        Coronagraph Inner Working Angle ( * lambda / diam)
    OWA : float, optional
        Coronagraph Outer Working Angle ( * lambda / diam)
    Tsys  : float, optional
        Telescope mirror temperature [K]
    Tdet  : float, optional
        Telescope detector temperature [K]
    emis : float, optional
        Effective emissivity for the observing system (of order unity)
    De : float, optional
        Dark current [counts/s]
    DNHpix : float, optional
        Number of horizontal/spatial pixels for dispersed spectrum
    Re : float, optional
        Read noise counts per pixel
    Dtmax : float, optional
        Detector maximum exposure time [hours]
    X : float, optional
        Width of photometric aperture ( * lambda / diam)
    qe : float, optional
        Detector quantum efficiency
    MzV : float, optional
        V-band zodiacal light surface brightness [mag/arcsec**2]
    MezV : float, optional
        V-band exozodiacal light surface brightness [mag/arcsec**2]
    wantsnr : float, optional
        Desired signal-to-noise ratio in each pixel
    FIX_OWA : bool, optional
        Set to fix OWA at OWA*lammin/D, as would occur if lenslet array is limiting the OWA
    COMPUTE_LAM : bool, optional
        Set to compute lo-res wavelength grid, otherwise the grid input as variable 'lam' is used
    SILENT : bool, optional
        Set to suppress print statements
    NIR : bool, optional
        Re-adjusts pixel size in NIR, as would occur if a second instrument was designed to handle the NIR
    THERMAL : bool, optional
        Set to compute thermal photon counts due to telescope temperature
    GROUND : bool, optional
        Set to simulate ground-based observations through atmosphere

    Returns
    -------
    lam : ndarray
        Observational wavelength grid [um]
    dlam : ndarray
        Observational spectral element width [um]
    A : ndarray
        Planetary geometric albedo on lam grid
    q : ndarray
        Quantum efficiency grid
    Cratio : ndarray
        Planet-star contrast ratio
    cp : ndarray
        Planetary photon count rate on detector
    csp : ndarray
        Speckle photon count rate on detector
    cz : ndarray
        Zodiacal photon count rate on detector
    cez : ndarray
        Exozodiacal photon count rate on detector
    cD : ndarray
        Dark current photon count rate on detector
    cR : ndarray
        Read noise photon count rate on detector
    cth : ndarray
        Instrument thermal photon count rate on detector
    DtSNR : ndarray
        Exposure time required to get desired S/N (wantsnr) [hours]
    """


    #check if user wants different Res_UV and Res_NIR:
    if Res_UV == -1: Res_UV = Res
    if Res_NIR == -1: Res_NIR = Res
        
    #Phi = 1.
    convolution_function = downbin_spec
    #convolution_function = degrade_spec

    # Configure for different telescope observing modes
    if mode == 'Imaging':
        filters = filter_wheel
        IMAGE = True
        COMPUTE_LAM = False
        # sorted filter dict by bandcenters
        tdict = sorted(filters.__dict__.iteritems(), key=lambda x: x[1].bandcenter)
        # Construct array of wavelengths
        lam = np.array([x[1].bandcenter for x in tdict])
        # Construct array of wavelength bin widths (FWHM)
        dlam = np.array([x[1].FWHM for x in tdict])
        Nlam = len(lam)
    elif mode == 'IFS':
        IMAGE = False
        COMPUTE_LAM = True
    else:
        print "Invalid telescope observing mode. Select 'IFS', or 'Imaging'."
        sys.exit()

    # fraction of planetary signal in Airy pattern
    fpa = f_airy(X)


    # Set wavelength grid
    if COMPUTE_LAM:
        lam, dlam = construct_lam(lammin, lammax, Res, UV=UV, NIR=NIR, Res_UV = Res_UV, Res_NIR = Res_NIR)
    elif IMAGE:
        pass
    else:
        # Throw error
        print "Error in make_noise: Not computing wavelength grid or providing filters!"
        return None

    # Set Quantum Efficiency
    q = set_quantum_efficiency(lam, qe, NIR=NIR)

    # Set Dark current and Read noise
    De = set_dark_current(lam, De_UV, De_VIS, De_NIR, lammax, Tdet, NIR=NIR, lammin_uv=lammin_uv, lammin_vis=lammin_vis, lammin_nir=lammin_nir)
    Re = set_read_noise(lam, Re_UV, Re_VIS, Re_NIR, NIR=NIR, lammin_uv=lammin_uv, lammin_vis=lammin_vis, lammin_nir=lammin_nir)

    # Set Angular size of lenslet
    theta = set_lenslet(lam, lammin, diam, NIR=NIR, UV=UV, lammin_vis = lammin_vis, lammin_nir=lammin_nir, lammin_uv=lammin_uv)

    # Set throughput
    #sep  = r/d*np.sin(alpha*np.pi/180.)*np.pi/180./3600. # separation in radians
   # T = set_throughput(lam, Tput, diam, sep, IWA, OWA, lammin, FIX_OWA=FIX_OWA, SILENT=SILENT)
    T = Tput
    # Modify throughput by atmospheric transmission if GROUND-based
    if GROUND:
        Tatmos = set_atmos_throughput(lam, dlam, convolution_function)
        # Multiply telescope throughput by atmospheric throughput
        T = T * Tatmos

    # Degrade albedo and stellar spectrum
    if COMPUTE_LAM:
        Fo = convolution_function(Fohr,lamhr,lam,dlam=dlam)
       # Fs = convolution_function(solhr, lamhr, lam, dlam=dlam)
    elif IMAGE:
        # Convolve with filter response
        Fo = convolve_spec(Fohr, lamhr, filters)
       # Fs = convolve_spec(solhr, lamhr, filters)
    else:
        Fo = Fohr
       # Fs = solhr

    # Compute fluxes
    #Fs = Fstar(lam, Teff, Rs, r, AU=True) # stellar flux on planet
    Fo = Fobj(Fo, Ro, d)         # planet flux at telescope
   # Cratio = FpFs(A, Phi, Rp, r)

    ##### Compute count rates #####
    co     =  cobj(q, fpa, T, lam, dlam, Fo, diam)
   # cp     =  cplan(q, fpa, T, lam, dlam, Fp, diam)                            # planet count rate
    cz     =  czodi(q, X, T, lam, dlam, diam, MzV)                           # solar system zodi count rate
   # cez    =  cezodi(q, X, T, lam, dlam, diam, r, \
      #  Fstar(lam,Teff,Rs,1.,AU=True), Nez, MezV)                            # exo-zodi count rate
   # csp    =  cspeck(q, T, C, lam, dlam, Fstar(lam,Teff,Rs,d), diam)         # speckle count rate
    cD     =  cdark(De, X, lam, diam, theta, DNHpix, IMAGE=IMAGE)            # dark current count rate
    cR     =  cread(Re, X, lam, diam, theta, DNHpix, Dtmax, IMAGE=IMAGE)     # readnoise count rate
    if THERMAL:
        cth    =  ctherm(q, X, lam, dlam, diam, Tsys, emis)                      # internal thermal count rate
    else:
        cth = np.zeros_like(cp)
    # Add earth thermal photons if GROUND
    if GROUND:
        # Compute ground intensity due to sky background
        Itherm  = get_thermal_ground_intensity(lam, dlam, convolution_function)
        # Compute Earth thermal photon count rate
        cthe = ctherm_earth(q, X, lam, dlam, diam, Itherm)
        # Add earth thermal photon counts to telescope thermal counts
        cth = cth + cthe
        if False:
            import matplotlib.pyplot as plt; from matplotlib import gridspec
            fig2 = plt.figure(figsize=(8,6))
            gs = gridspec.GridSpec(1,1)
            ax1 = plt.subplot(gs[0])
            ax1.plot(lam, cthe, c="blue", ls="steps-mid", label="Earth Thermal")
            ax1.plot(lam, cth, c="red", ls="steps-mid", label="Telescope Thermal")
            ax1.plot(lam, cp, c="k", ls="steps-mid", label="Planet")
            ax1.set_ylabel("Photon Count Rate [1/s]")
            ax1.set_xlabel("Wavelength [um]")
            plt.show()

    cb = (cz + cD + cR + cth)
    cnoise =  co + 2*cb                # assumes background subtraction
    ctot = co + cz + cD + cR + cth

    '''
    Giada: where does the factor of 2 come from [above]?

    Ty (Via email): That's due to "background subtraction".
    If you were to take a single exposure, and had the ability
    to post-process the image to the Poisson noise limit, you
    wouldn't have the factor of two.  However, it's not yet
    clear that we'll be able to reach the Poisson, single observation limit.
    Instead, the current idea is that you take two observations (with
    exposure time Delta t/2), with the telescope rotated by a
    small amount between exposures, and then subtract the two images.
    So, for a fixed exoplanet count (cp), the roll technique would
    give you 2x as many noise counts due to background sources as
    would the single-observation technique.
    See also the top of page 4 of Brown (2005).
    '''

    # Exposure time to SNR
    DtSNR = exptime_element(lam, co, cnoise, wantsnr)
    print Fo
    return lam, dlam, Fo, q, co, cz, cD, cR, cth, DtSNR
