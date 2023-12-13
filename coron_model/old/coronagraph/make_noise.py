# Import dependent modules
import numpy as np
from .degrade_spec import degrade_spec, downbin_spec
from .convolve_spec import convolve_spec
from .noise_routines import Fstar, Fplan, FpFs, cplan, czodi, cezodi, cspeck, cdark, cread, ctherm, ccic, f_airy
import pdb

__all__ = ['make_noise']

def make_noise(Ahr, lamhr, solhr, telescope, planet, star, wantsnr=10.0, FIX_OWA = False, COMPUTE_LAM = False,\
    SILENT = False, NIR = True, THERMAL = False):
    """
    Generate photon count rates for specified telescope and planet parameters

    Parameters
    ----------
    Ahr : array
        hi-res planetary albedo spectrum
    lamhr : array
        wavelength grid for Ahr (um)
    solhr : array
        hi-res TOA solar spectrum (W/m**2/um)
    telescope : Telescope
        Telescope object containing parameters
    planet : Planet
        Planet object containing parameters
    star : Star
        Star object containing parameters
    FIX_OWA : bool
        set to fix OWA at OWA*lammin/D, as would occur if lenslet array is limiting the OWA
    COMPUTE_LAM : bool
        set to compute lo-res wavelength grid, otherwise the grid input as variable 'lam' is used
    NIR : bool
        re-adjusts pixel size in NIR, as would occur if a second instrument was designed to handle the NIR
    THERMAL : bool
        set to compute thermal photon counts due to telescope temperature
    """

    # Assign Class attributes to local variables
    alpha = planet.alpha
    Phi = planet.Phi
    Rp = planet.Rp
    Teff = star.Teff
    Rs=star.Rs
    r=planet.a
    d=planet.distance
    Nez=planet.Nez
    lammin=telescope.lammin
    lammax=telescope.lammax
    Res=telescope.resolution
    diam=telescope.diameter
    Tput=telescope.throughput
    C=telescope.contrast
    IWA=telescope.IWA
    OWA=telescope.OWA
    Tsys=telescope.Tsys
    Tdet=telescope.Tdet
    emis=telescope.emissivity

    # Set the Imaging Mode?
    if telescope.mode == 'Imaging':
        filters = telescope.filter_wheel
        IMAGE = True
        COMPUTE_LAM = False
        # sorted filter dict by bandcenters
        tdict = sorted(filters.__dict__.iteritems(), key=lambda x: x[1].bandcenter)
        # Construct array of wavelengths
        lam = np.array([x[1].bandcenter for x in tdict])
        # Construct array of wavelength bin widths (FWHM)
        dlam = np.array([x[1].FWHM for x in tdict])
        Nlam = len(lam)
    else:
        IMAGE = False

    # Set key system parameters
    De     = telescope.darkcurrent  # dark current (s**-1)
    DNHpix = telescope.DNHpix       # horizontal pixel spread of IFS spectrum
    Re     = telescope.readnoise    # read noise per pixel
    Dtmax  = telescope.Dtmax        # maximum exposure time (hr)
    X      = telescope.X            # size of photometric aperture (lambda/D)
    qe     = telescope.qe           # quantum efficiency

    # fraction of planetary signal in Airy pattern
    fpa = f_airy(X)

    # Set astrophysical parameters
    MzV  = 23.0     # zodiacal light surface brightness (mag/arcsec**2)
    MezV = 22.0     # exozodiacal light surface brightness (mag/arcsec**2)

    # Set wavelength grid
    if COMPUTE_LAM:
        lam  = lammin #in [um]
        Nlam = 1
        while (lam < lammax):
            lam  = lam + lam/Res
            Nlam = Nlam +1
        lam    = np.zeros(Nlam)
        lam[0] = lammin
        for j in range(1,Nlam):
            lam[j] = lam[j-1] + lam[j-1]/Res
        Nlam = len(lam)
        dlam = np.zeros(Nlam) #grid widths (um)
        # Set wavelength widths
        for j in range(1,Nlam-1):
            dlam[j] = 0.5*(lam[j+1]+lam[j]) - 0.5*(lam[j-1]+lam[j])
        #widths at edges are same as neighbor
        dlam[0] = dlam[1]
        dlam[Nlam-1] = dlam[Nlam-2]
    elif IMAGE:
        pass
    else:
        # Throw error
        print "Error in make_noise: Not computing wavelength grid or providing filters!"
        return None

    # Set Quantum Efficiency
    q = np.zeros(Nlam) + qe
   # for j in range(Nlam):
   #     if (lam[j] <= 0.7):
   #         q[j] = qe
   #     else:
   #         q[j] = qe*(1.0 - (lam[j]-0.7)/(1.0-0.7))
   # if q[j] < 0.0: q[j] = 0.

    # Set Dark current and Read noise
    De = np.zeros(Nlam) + De
    Re = np.zeros(Nlam) + Re

    # Set Angular size of lenslet
    theta = lammin/1.e6/diam/2.*(180/np.pi*3600.) #assumes sampled at ~lambda/2D (arcsec)
    if NIR:
        theta = np.zeros(Nlam)
        iVIS  = (lam <= 1.0)
        iNIR  = (lam > 1.0)
        theta[iVIS] = lammin/1e6/diam/2.*(180/np.pi*3600.)
        theta[iNIR] = 1.0/1e6/diam/2.*(180/np.pi*3600.)
      #  q[iNIR]  = qe
        Re[iNIR] = 2.                # Different from vis detector
        De[iNIR] = 1e-3              # Different from vis detector

    # Set throughput
    T    = Tput + np.zeros(Nlam)
    sep  = r/d*np.sin(alpha*np.pi/180.)*np.pi/180./3600. # separation in radians
    iIWA = ( sep < IWA*lam/diam/1.e6 )
    if (True if True in iIWA else False):
        T[iIWA] = 0. #zero transmission for points inside IWA have no throughput
        if ~SILENT:
            print 'WARNING: portions of spectrum inside IWA'
    if FIX_OWA:
        if ( sep > OWA*lammin/diam/1.e6 ):
            T[:] = 0. #planet outside OWA, where there is no throughput
            if ~SILENT:
                print 'WARNING: planet outside fixed OWA'
    else:
        iOWA = ( sep > OWA*lam/diam/1.e6 )
        if (True if True in iOWA else False):
            T[iOWA] = 0. #points outside OWA have no throughput
            if ~SILENT:
                print 'WARNING: portions of spectrum outside OWA'


    # Degrade albedo and stellar spectrum
    if COMPUTE_LAM:
        A = degrade_spec(Ahr,lamhr,lam,dlam=dlam)
        Fs = degrade_spec(solhr, lamhr, lam, dlam=dlam)
    elif IMAGE:
        # Convolve with filter response
        A = convolve_spec(Ahr, lamhr, filters)
        Fs = convolve_spec(solhr, lamhr, filters)
    else:
        A = Ahr
        Fs = solhr

    # Compute fluxes
    #Fs = Fstar(lam, Teff, Rs, r, AU=True) # stellar flux on planet
    Fp = Fplan(A, Phi, Fs, Rp, d)         # planet flux at telescope
    Cratio = FpFs(A, Phi, Rp, r)


    ##### Compute count rates #####
    cp     =  cplan(q, fpa, T, lam, dlam, Fp, diam)                            # planet count rate
    cz     =  czodi(q, X, T, lam, dlam, diam, MzV)                           # solar system zodi count rate
    cez    =  cezodi(q, X, T, lam, dlam, diam, r, \
        Fstar(lam,Teff,Rs,1.,AU=True), Nez, MezV)                            # exo-zodi count rate
    csp    =  cspeck(q, T, C, lam, dlam, Fstar(lam,Teff,Rs,d), diam)         # speckle count rate
    cD     =  cdark(De, X, lam, diam, theta, DNHpix, IMAGE=IMAGE)            # dark current count rate
    cR     =  cread(Re, X, lam, diam, theta, DNHpix, Dtmax, IMAGE=IMAGE)     # readnoise count rate
    if THERMAL:
        cth    =  ctherm(q, X, lam, dlam, diam, Tsys, emis)                      # internal thermal count rate
    else:
        cth = np.zeros_like(cp)
    cnoise =  cp + 2*(cz + cez + csp + cD + cR + cth)                        # assumes background subtraction
    cb = (cz + cez + csp + cD + cR + cth)
    ctot = cp + cz + cez + csp + cD + cR + cth

    '''
    Giada: where does the factor of 2 come from?

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
    DtSNR = np.zeros(Nlam)
    DtSNR[:] = 0.
    i = (cp > 0.)
    if (True if True in i else False):
        DtSNR[i] = (wantsnr**2.*cnoise[i])/cp[i]**2./3600. # (hr)
    # added by Giada:
    #if whichplanet == 'earth':
    #    print 'Functionality not added to python version... yet.'
        #pt5 = closest(lam, 0.55) ;same wavelength chris stark used
        #time = dtsnr(pt5)*3600.*1
        #save, time, filename='~/idl/noise_model/earthtime.sav'
    #if whichplanet != 'earth':
    #    print 'Functionality not added to python version... yet.'
        #then restore, '~/idl/noise_model/earthtime.sav'

    # These pieces are fundamental, but should go outside this function
    # as they depend on the particular exposure time, not just the telescope
    #noisyspec = np.random.poisson(cnoise * time)
    #planet = noisyspec - 2.*(cz + cez + csp + cD + cR + cth)*  time
    #sun = (time * cp)/A
    #SNR = cp*time/np.sqrt((cp+2*cb)*time)
    #noisy = np.random.randn(len(A))*A/SNR+A

    return lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, DtSNR
