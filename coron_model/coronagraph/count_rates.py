# Import dependent modules
import numpy as np
import sys
from .degrade_spec import downbin_spec
from .convolve_spec import convolve_spec
from .noise_routines import Fstar, Fplan, FpFs, cplan, czodi, cezodi, cspeck, \
    cdark, cread, ctherm, f_airy, ctherm_earth, construct_lam, \
    set_quantum_efficiency, set_read_noise, set_dark_current, set_cic, set_lenslet, \
    set_throughput, set_atmos_throughput, get_thermal_ground_intensity, \
    exptime_element, lambertPhaseFunction


__all__ = ['count_rates']

def count_rates(Ahr, lamhr, solhr,
                alpha, Rp, Teff, Rs, r, d, Nez,
                mode   = "IFS",
                filter_wheel = None,
                Phi = -1, #chosen as a value a user would never set 
                lammin = 0.55,
                lammax = 2.,
                Res    = 140.0,
                diam   = 15.,
                diam_ins = -1,
                collect_area = -1,
                Tput   = 0.15,
                Tput_uv = 0.15,
                Tput_nir = 0.15,
                o_Tput_uv = 0.12,
                o_Tput_vis = 0.32,
                o_Tput_nir = 0.60,
                C      = 1e-10,
                IWA    = 3.5,
                IWA_UV = -1,
                IWA_NIR = -1,
                OWA    = 64.0,
                OWA_UV = -1,
                OWA_NIR = -1,
                Tsys   = 150.0,
                Tdet   = 50.0,
                emis   = 0.09,
                De_UV  = 3e-5, #BOL dark current 
                De_VIS = 3e-5, #BOL dark current
                De_NIR = 2e-3, #specs for teledyne detector
                CIC_UV = 1.3e-3,
                CIC_VIS = 1.3e-3,
                CIC_NIR = 0.0,
                DNHpix = 3.0,
                Re_UV  = 0.0,
                Re_VIS = 0.0,
                Re_NIR = 2.5,
                Dtmax  = 1.0,
                X      = 0.7, #changed from 1.5 
                qe     = 0.675, #QE (0.9) times charge transfer term (0.75)
                MzV    = 23.0, #23
                MezV   = 22.0, #22
                Res_NIR = -1, 
                Res_UV = -1,
                lammin_uv = 0.2,
                lammin_vis = 0.515,
                lammin_nir = 1.0,
                ntherm = 1,
                gain = 1,
                bg_factor = 2., #2 assumes roll manouver for noise calculation
                mirror = 'perfect', #choices = perfect, Al, Au
                ssIWA = -1, #switched off if negative
                ssOWA = -1, #switched off if negative
                wantsnr=10.0, FIX_OWA = False, COMPUTE_LAM = False,
                SILENT = False, NIR = True, UV=True, THERMAL = True,
                GROUND = False, writeout=False, writeoutpath='',
                rebin_vis=False, LUVOIR_A = False):
    """
    Runs coronagraph model (Robinson et al., 2016) to calculate planet and noise
    photon count rates for specified telescope and system parameters.

in    Parameters
    ----------
    Ahr : array
        High-res, wavelength-dependent planetary geometric albedo
    lamhr : arrays
        High-res wavelength grid  [um]
    solhr : array
        High-res TOA solar spectrum [W/m**2/um]
    alpha : float
        Planet phase angle [deg]
    Phi : float
        Planet phase function [deg] (will default to calculating with alpha if not set)
    Rp : float
        Planet radius [R_earth]
    Teff : float
        Stellar effective temperature [K]
    Rs : float
        Stellar radius [R_sun]
    r : float
        Planet semi-major axis [AU]
    d : float
        Distance to observed star-planet system [pc]
    Nez : float
        Number of exozodis in exoplanetary disk
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
    collect area : float, optional
        If user sets this, they can specify a collecting area different from the diameter [m^2]
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
    ntherm : float, optional
        Number of thermal surfaces
    ssIWA : float, optional
        Starshade-like IWA (fixed angle in mas)
    ssOWA : float, optional
        Starshade-like OWA (fixed angle in arcsec)
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
    if Phi == -1: #check if user has manually set Phi (they'd never set it to -1)
        #calculate phase function
        Phi = lambertPhaseFunction(alpha)

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
        # Nlam = len(lam)
    elif mode == 'IFS':
        IMAGE = False
        COMPUTE_LAM = True
    else:
        print("Invalid telescope observing mode. Select 'IFS', or 'Imaging'.") 
        sys.exit()

    # fraction of planetary signal in Airy pattern
    fpa = f_airy(X)

    
    #check if user has specified a different collecting area
    if collect_area == -1:
        diam2 = diam
    if collect_area != -1:
        diam2 = np.sqrt(collect_area/3.14159)*2.

    #check if user has specified an inscribed circle
    if diam_ins == -1:
        diam_ins = diam
 
        
    # Set wavelength grid
    if COMPUTE_LAM:
        lam, dlam = construct_lam(lammin, lammax, Res, UV=UV, NIR=NIR, Res_UV = Res_UV, Res_NIR = Res_NIR, lammin_uv=lammin_uv, lammin_vis=lammin_vis, lammin_nir=lammin_nir)
    elif IMAGE:
        pass
    else:
        # Throw error
        print("Error in make_noise: Not computing wavelength grid or providing filters!") 
        return None

    # Set Quantum Efficiency
    q = set_quantum_efficiency(lam, qe, NIR=NIR)

    # Set Dark current and Read noise and clock induced charge
    print(De_UV, De_VIS, De_NIR) 
    De = set_dark_current(lam, De_UV, De_VIS, De_NIR, lammax, Tdet, NIR=NIR, lammin_uv=lammin_uv, lammin_vis=lammin_vis, lammin_nir=lammin_nir)
    Re = set_read_noise(lam, Re_UV, Re_VIS, Re_NIR, NIR=NIR, lammin_uv=lammin_uv, lammin_vis=lammin_vis, lammin_nir=lammin_nir)
    cic = set_cic(lam, CIC_UV, CIC_VIS, CIC_NIR, NIR=NIR, lammin_uv=lammin_uv, lammin_vis=lammin_vis, lammin_nir=lammin_nir)
    
    # Set Angular size of lenslet
    theta = set_lenslet(lam, lammin, diam_ins, NIR=NIR, UV=UV, lammin_vis = lammin_vis, lammin_nir=lammin_nir, lammin_uv=lammin_uv)
    
    # Set throughput
    sep  = r/d*np.sin(alpha*np.pi/180.)*np.pi/180./3600. # separation in radians


    if ssIWA != -1:
        ssIWArad = ssIWA * (np.pi/648000000.) # ssIWA in radians
    else:
        ssIWArad = -1
    if ssIWArad != -1:
        ssOWArad = ssIWA * (np.pi/648000.) # ssIWA in radians
    else:
        ssOWArad = -1
    T = set_throughput(lam, Tput, Tput_uv, Tput_nir, o_Tput_vis, o_Tput_uv, o_Tput_nir, diam, sep, IWA, IWA_UV, IWA_NIR, OWA, OWA_UV, OWA_NIR, ssIWArad, ssOWArad, lammin, mirror, ntherm,FIX_OWA=FIX_OWA, SILENT=SILENT, lammin_uv=lammin_uv, lammin_vis=lammin_vis, lammin_nir=lammin_nir, LUVOIR_A = LUVOIR_A)


    
    # Modify throughput by atmospheric transmission if GROUND-based
    if GROUND:
        Tatmos = set_atmos_throughput(lam, dlam, convolution_function)
        # Multiply telescope throughput by atmospheric throughput
        T = T * Tatmos

    # Degrade albedo and stellar spectrum
    if COMPUTE_LAM:
        A = convolution_function(Ahr,lamhr,lam,dlam=dlam)
        Fs = convolution_function(solhr, lamhr, lam, dlam=dlam)
    elif IMAGE:
        # Convolve with filter response
        A = convolve_spec(Ahr, lamhr, filters)
        Fs = convolve_spec(solhr, lamhr, filters)
    else:
        A = Ahr
        Fs = solhr

    # Compute fluxes
    Fs = Fstar(lam, Teff, Rs, r, AU=True) # stellar flux on planet
    Fp = Fplan(A, Phi, Fs, Rp, d)         # planet flux at telescope
    # i =  (lam >= 0.54) & (lam <= 0.56)

    
    Cratio = FpFs(A, Phi, Rp, r)

    ##### Compute count rates #####
    cp     =  cplan(q, fpa, T, lam, dlam, Fp, diam2)                            # planet count rate
    cp = cp * gain
    cz     =  czodi(q, X, T, lam, dlam, diam2, MzV)                           # solar system zodi count rate
    cz = cz * gain
    cez    =  cezodi(q, X, T, lam, dlam, diam2, r, \
        Fstar(lam,Teff,Rs,1.,AU=True), Nez, MezV)                            # exo-zodi count rate
    cez = cez * gain
    csp    =  cspeck(q, T, C, lam, dlam, Fstar(lam,Teff,Rs,d), diam2)         # speckle count rate
    csp = csp * gain
    cD     =  cdark(De, X, lam, diam2, theta, DNHpix, IMAGE=IMAGE)            # dark current count rate
    cR     =  cread(Re, X, lam, diam2, theta, DNHpix, Dtmax, IMAGE=IMAGE)     # readnoise count rate

    if THERMAL:
        cth    =  ntherm*ctherm(q, X, lam, dlam, diam2, Tsys, emis, T)                      # internal thermal count rate
        cth = cth * gain
    else:
        cth = np.zeros_like(cp)
    # Add earth thermal photons if GROUND
    if GROUND:
        # Compute ground intensity due to sky background
        Itherm  = get_thermal_ground_intensity(lam, dlam, convolution_function)
        # Compute Earth thermal photon count rate
        cthe = ctherm_earth(q, X, lam, dlam, diam2, Itherm)
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
            
    ctot = cp + cz + cez + csp + cD + cR + cth
    ccic = ctot * cic #is this right??
    ctot = ctot + ccic

    cb = (cz + cez + csp + cD + cR + cth + ccic)
    cnoise =  cp + bg_factor*cb                # assumes background subtraction
    ctot = cp + cz + cez + csp + cD + cR + cth
   # ccic = ctot * cic #is this right??
   # ctot = ctot + ccic
   # cb = cb + ccic 

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
    DtSNR = exptime_element(lam, cp, cnoise, wantsnr)

    #if rebin_vis then rebin visible channel to R=75 for everywhere except O2 A band
    if rebin_vis: 
        Res_lo = 75
        lam_lo, dlam_lo = construct_lam(lammin, lammax, Res_lo, UV=UV, NIR=NIR, Res_UV = Res_UV, Res_NIR = Res_NIR, lammin_uv=lammin_uv, lammin_vis=lammin_vis, lammin_nir=lammin_nir)
        #stitch wl grids together around oxygen:
        io2 = (lam > 0.73) & (lam < 0.8)
        ishort = (lam_lo <= 0.73)
        ilong = (lam_lo >= 0.8)
        lam_rebin = np.concatenate([lam_lo[ishort], lam[io2], lam_lo[ilong]])
        dlam_rebin = np.concatenate([dlam_lo[ishort], dlam[io2], dlam_lo[ilong]])
        #downbin everything now 
        cp_lo = convolution_function(cp,lam,lam_rebin,dlam=dlam_rebin)
        cnoise_lo = convolution_function(cnoise,lam,lam_rebin,dlam=dlam_rebin)
        DtSNR_lo = exptime_element(lam_rebin, cp_lo, cnoise_lo, wantsnr)

        

    
    
    # write text output file
    if writeout:
        data_tag = writeoutpath+'output.txt'
        y_sav = np.array([lam,Cratio,A,q,cp,csp,cz,cez,cD,cR,cth,DtSNR])
        np.savetxt(data_tag, y_sav.T)
        print('Saved: ' + data_tag) 

        print(lam) 
    return lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, ccic, DtSNR
