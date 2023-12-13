import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import rc
import os
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rcParams['font.size'] = 20.0

from .make_noise import make_noise
from .teleplanstar import Telescope, Planet, Star

__all__ = ['generate_observation', 'smart_observation', 'planetzoo_observation', 'process_noise', 'exptime_band', 'interp_cont_over_band']

planetdir = "planets/"
relpath = os.path.join(os.path.dirname(__file__), planetdir)

def planetzoo_observation(name='earth', telescope=Telescope(), planet=Planet(), itime=10.0,
                            planetdir = relpath, plot=True, savedata=False, saveplot=False,
                            ref_lam=0.55, THERMAL=False):
    """Uses coronagraph model to observe planets located in planetdir

    Parameters
    ----------
    name : str (optional)
        Name of the planet
    telescope : Telescope (optional)
        Telescope object to be used for observation
    planet : Planet (optional)
        Planet object to be used for observation
    itime : float (optional)
        Integration time (hours)
    planetdir : str
        Location of planets/ directory
    plot : bool (optional)
        Make plot flag
    savedata : bool (optional)
        Save output as data file
    saveplot : bool (optional)
        Save plot as PDF
    ref_lam : float (optional)
        Wavelength at which SNR is computed

    Returns
    -------
    lam : array
        Observed wavelength array (microns)
    spec : array
        Observed reflectivity spectrum
    sig : array
        Observed 1-sigma error bars on spectrum
    """

    '''
    planet choices:
        earth, venus, archean,
        earlymars, hazyarchean, earlyvenus,
        jupiter, saturn, uranus, neptune, mars,
        fstarozone
    star choices (for now):
    sun, f2v
    '''

    import os
    try:
        l = os.listdir(planetdir)
    except OSError:
        print "Error in planetzoo_observation(): planetdir does not exist in current location. \nSet planetdir='location/of/planets/'"
        return None

    whichplanet = name
    startype = planet.star
    tag = name

    if True:

        #if whichplanet == 'earth':
        fn = 'earth_quadrature_radiance_refl.dat'
        model = np.loadtxt(os.path.join(planetdir,fn), skiprows=8)
        lamhr = model[:,0]
        radhr = model[:,1]
        solhr = model[:,2]
        reflhr = model[:,3]
        Ahr   = np.pi*(np.pi*radhr/solhr) # hi-resolution reflectivity
        planet.Rp    = 1.0     # Earth radii
        planet.r     = 1.0     # semi-major axis (AU)

        if whichplanet == 'venus':
            fn = 'Venus_geo_albedo.txt'
            model = np.loadtxt(os.path.join(planetdir,fn))
            lamhr = model[:,0]
            Ahr   = model[:,1]
            planet.Rp    = 0.95     #Earth radii
            planet.r     = 0.72     #semi-major axis (AU)

        if whichplanet == 'archean':
            fn = 'ArcheanEarth_geo_albedo.txt'
            model = np.loadtxt(os.path.join(planetdir,fn))
            lamhr = model[:,0]
            Ahr   = model[:,1]
            planet.Rp    = 1.0     #Earth radii
            planet.r     = 1.0     #semi-major axis (AU)

        if whichplanet == 'earlymars':
            fn = 'EarlyMars_geo_albedo.txt'
            model = np.loadtxt(os.path.join(planetdir,fn))
            lamhr = model[:,0]
            Ahr   = model[:,1]
            planet.Rp    = 0.53     #Earth radii
            planet.r     = 1.52     #semi-major axis (AU)

        if whichplanet == 'hazyarchean':
            fn = 'Hazy_ArcheanEarth_geo_albedo.txt'
            model = np.loadtxt(os.path.join(planetdir,fn))
            lamhr = model[:,0]
            Ahr   = model[:,1]
            planet.Rp    = 1.0     #Earth radii
            planet.r     = 1.0     #semi-major axis (AU)

        if whichplanet == 'earlyvenus':
            fn = 'EarlyVenus_geo_albedo.txt'
            model = np.loadtxt(os.path.join(planetdir,fn))
            lamhr = model[:,0]
            Ahr   = model[:,1]
            planet.Rp    = 0.95     #Earth radii
            planet.r     = 0.72     #semi-major axis (AU)

        if whichplanet == 'jupiter':
            fn = planetdir+'Jupiter_geo_albedo.txt'
            model = np.loadtxt(os.path.join(planetdir,fn))
            lamhr = model[:,0]
            Ahr   = model[:,1]
            planet.Rp    = 10.86     #Earth radii
            planet.r     = 5.20     #semi-major axis (AU)

        if whichplanet == 'saturn':
            fn = planetdir+'Saturn_geo_albedo.txt'
            model = np.loadtxt(os.path.join(planetdir,fn))
            lamhr = model[:,0]
            Ahr   = model[:,1]
            planet.Rp    = 9.00     #Earth radii
            planet.r     = 9.54     #semi-major axis (AU)

        if whichplanet == 'uranus':
            fn = 'Uranus_geo_albedo.txt'
            model = np.loadtxt(os.path.join(planetdir,fn))
            lamhr = model[:,0]
            Ahr   = model[:,1]
            planet.Rp    = 3.97     #Earth radii
            planet.r     = 19.19     #semi-major axis (AU)


        if whichplanet == 'warmuranus':
            fn = 'Uranus_geo_albedo.txt'
            model = np.loadtxt(os.path.join(planetdir,fn))
            lamhr = model[:,0]
            Ahr   = model[:,1]
            planet.Rp    = 3.97     #Earth radii
            planet.r     = 5.20     #semi-major axis (AU)

        if whichplanet == 'warmneptune':
            fn = 'Neptune_geo_albedo.txt'
            model = np.loadtxt(os.path.join(planetdir,fn))
            lamhr = model[:,0]
            Ahr   = model[:,1]
            planet.Rp    = 3.97     #Earth radii
            planet.r     = 5.20     #semi-major axis (AU)

        if whichplanet == 'neptune':
            fn = 'Neptune_geo_albedo.txt'
            model = np.loadtxt(os.path.join(planetdir,fn))
            lamhr = model[:,0]
            Ahr   = model[:,1]
            planet.Rp    = 3.85     #Earth radii
            planet.r     = 30.07     #semi-major axis (AU)


        if whichplanet == 'mars':
            fn = 'Mars_geo_albedo.txt'
            model = np.loadtxt(os.path.join(planetdir,fn))
            lamhr = model[:,0]
            Ahr   = model[:,1]
            planet.Rp    = 0.53     #Earth radii
            planet.r     = 1.52     #semi-major axis (AU)

        if whichplanet == 'fstarozone':
            print('fstarozone functionality not yet added')
            #fn = 'fstarcloudy.sav'
            #fn = 'F2V_5.e-1fCO2_1.e6H2Volc_1.e10BIF.out_toa.rad'
            #READCOL, 'planets/'+fn, lamhr, wnhr, sunhr, fluxhr
            #restore, 'planets/'+fn
            #lamhr = reverse(wl_)
            #ahr = reverse(cloudyalb_)
            #Ahr = (2/3.) * fluxhr/(sunhr/2.) #geometric albedo & phase corrections
            #Rp    = 1.0     #Earth radii
            #r     = 1.72    #semi-major axis (AU)



    # star parameters
    if startype == '':
        Teff  = 5780.   #Sun-like Teff (K)
        Rs    = 1.      #star radius in solar radii
    if  startype == 'sun':
        Teff  = 5780.   #Sun-like Teff (K)
        Rs    = 1.      #star radius in solar radii
    if  startype == 'f2v':
        Teff  = 7050.   #Sun-like Teff (K)
        Rs    = 1.3      #star radius in solar radii

    star = Star(Teff=Teff, Rs=Rs)

    # Shawn: "I don't like noise.  It makes me sad."

    lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, DtSNR \
        = make_noise(Ahr, lamhr, solhr, telescope, planet, star, COMPUTE_LAM=True, THERMAL=THERMAL)

    # Calculate background photon count rate
    cb = (cz + cez + csp + cD + cR + cth)

    # Calculate the SNR of observation
    time = itime * 3600. # Convert hours to seconds
    #SNR = calc_SNR(time, cp, cb)

    # Generate noisy spectrum by drawing data points from a normal distribution
    #spec, sig = draw_noisy_spec(A, SNR)

    # Calculate SNR, sigma, and noised-up spectrum
    spec, sig, SNR = process_noise(time, Cratio, cp, cb)

    if plot:
        plot_coronagraph_spectrum(lam, spec, sig, itime, planet.distance, ref_lam, SNR, truth=Cratio)

    # Save Synthetic data file (wavelength, albedo, error) if savedata=True
    if savedata:
        data_tag = 'observed_'+tag+'.txt'
        y_sav = np.array([lam,spec,sig])
        np.savetxt(data_tag, y_sav.T)
        print 'Saved: '+data_tag

    # Return Synthetic data and high-res spec
    return lam, spec, sig

def generate_observation(wlhr, Ahr, solhr, itime, telescope, planet, star,
                         ref_lam=0.55, tag='', plot=True, saveplot=False, savedata=False,
                         THERMAL=False, wantsnr=10):
    """
    Parameters
    ----------
    wlhr : float
        Wavelength array (microns)
    Ahr : float
        Geometric albedo spectrum array
    itime : float
        Integration time (hours)
    telescope : Telescope
        Telescope object
    planet : Planet
        Planet object
    star : Star
        Star object
    tag : string
        ID for output files
    plot : boolean
        Set to True to make plot
    saveplot : boolean
        Set to True to save the plot as a PDF
    savedata : boolean
        Set to True to save data file of observation

    Returns
    -------
    lam : array
        Wavelength grid of observed spectrum
    dlam: array

    A : array

    spec : array
        Albedo grid of observed spectrum
    sig : array
        One sigma errorbars on albedo spectrum
    SNR : array


    Output
    ------
    If saveplot=True then plot will be saved
    If savedata=True then data will be saved
    """

    # Skip call_noise and just call: noise
    lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, DtSNR = \
        make_noise(Ahr, wlhr, solhr, telescope, planet, star, wantsnr=wantsnr, COMPUTE_LAM=True, THERMAL=THERMAL)

    # Calculate background photon count rate
    cb = (cz + cez + csp + cD + cR + cth)

    # Calculate the SNR of observation
    time = itime * 3600. # Convert hours to seconds

    # Calculate SNR, sigma, and noised-up spectrum
    spec, sig, SNR = process_noise(time, Cratio, cp, cb)

    #SNR = calc_SNR(time, cp, cb)

    # Generate noisy spectrum by drawing data points from a normal distribution
    #spec, sig = draw_noisy_spec(A, SNR)


    if plot:
        plot_coronagraph_spectrum(lam, spec, sig, itime, planet.distance, ref_lam, SNR, truth=Cratio)

    # Save Synthetic data file (wavelength, albedo, error) if savedata=True
    if savedata:
        data_tag = 'observed_'+tag+'.txt'
        y_sav = np.array([lam,spec,sig])
        np.savetxt(data_tag, y_sav.T)
        print 'Saved: '+data_tag

    # Return Synthetic data and high-res spec

    return lam, dlam, Cratio, spec, sig, SNR


def smart_observation(radfile, itime, telescope, planet, star,
                         ref_lam=0.55, tag='', plot=True, saveplot=False, savedata=False,
                         THERMAL=False, wantsnr=10.):
    """Uses coronagraph noise model to create an observation of high resolution SMART output.

    Parameters
    ----------
    radfile : string
        Location and name of file to be read in
    itime : float
        Integration time (hours)
    telescope : Telescope
        Telescope object
    planet : Planet
        Planet object
    star : Star
        Star object
    tag : string
        ID for output files
    plot : boolean
        Set to True to make plot
    saveplot : boolean
        Set to True to save the plot as a PDF
    savedata : boolean
        Set to True to save data file of observation

    Returns
    ----------
    lam : array
        Wavelength grid of observed spectrum
    spec : array
        Albedo grid of observed spectrum
    sig : array
        One sigma errorbars on albedo spectrum
    rwl : array
        Wavelength grid of SMART output
    Ahr : array
        Albedo grid of SMART output

    Output
    ---------
    If saveplot=True then plot will be saved
    If savedata=True then data will be saved
    """

    # try importing readsmart
    try:
        import readsmart as rs
    except ImportError:
        print "Module 'readsmart' not found. Please install on your local machine \
        to proceed with this function. The source can be found at: \
        https://github.com/jlustigy/readsmart"
        return None, None, None, None, None

    # Read-in .rad file
    wlhr, wno, solar_spec, TOA_flux, rad_streams = rs.rad(radfile,getdata=True)

    # Calculate Hi-res reflectivity spectrum
    Ahr = (TOA_flux / solar_spec) #* np.pi / planet.Phi

    # Possibly convolve with gaussian?

    # Skip call_noise and just call: noise
    lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, DtSNR = \
        make_noise(Ahr, wlhr, solar_spec, telescope, planet, star, wantsnr=wantsnr, COMPUTE_LAM=True, THERMAL=THERMAL)

    # Calculate background photon count rate
    cb = (cz + cez + csp + cD + cR + cth)

    # Calculate the SNR of observation
    time = itime * 3600. # Convert hours to seconds
    #SNR = calc_SNR(time, cp, cb)

    # Generate noisy spectrum by drawing data points from a normal distribution
    #spec, sig = draw_noisy_spec(A, SNR)

    # Calculate SNR, sigma, and noised-up spectrum
    spec, sig, SNR = process_noise(time, Cratio, cp, cb)

    if plot:
        plot_coronagraph_spectrum(lam, spec, sig, itime, planet.distance, ref_lam, SNR, truth=Cratio)

    # Save Synthetic data file (wavelength, albedo, error) if savedata=True
    if savedata:
        data_tag = 'observed_smart_'+tag+'.txt'
        y_sav = np.array([lam,spec,sig])
        np.savetxt(data_tag, y_sav.T)
        print 'Saved: '+data_tag

    # Return Synthetic data and high-res spec

    return lam, spec, sig, wlhr, Ahr

def plot_spec():
    """
    Plots noised-up spectrum.
    """

    # Set string for plot text
    if itime > 2.0:
        timestr = "{:.0f}".format(itime)+' hours'
    else:
        timestr = "{:.0f}".format(itime*60)+' mins'
    plot_text = r'Distance = '+"{:.1f}".format(planet.distance)+' pc'+\
    '\n Integration time = '+timestr

    # If a reference wavelength is specified then return the SNR at that wl
    # corresponding to the integration time given
    if ref_lam:
        ireflam = find_nearest(lam,ref_lam)
        ref_SNR = SNR[ireflam]
        plot_text = plot_text + '\n SNR = '+"{:.1f}".format(ref_SNR)+\
            ' at '+"{:.2f}".format(lam[ireflam])+r' $\mu$m'

    # Plot observed spectrum; save pdf if saveplot=True
    lammin,lammax = np.min(lam)-0.1, np.max(lam)+0.1
    Amin, Amax = np.min(A)-np.max(sig)*1.1, np.max(A)+np.max(sig)*1.1
    tmin = np.min(Ahr[(wlhr > lammin) & (wlhr < lammax)])
    tmax = np.max(Ahr[(wlhr > lammin) & (wlhr < lammax)])
    if tmin < Amin: Amin = tmin
    if tmax > Amax: Amax = tmax
    #ymin,ymax = np.min(A), np.max(A)
    plot_tag = 'observed_smart_'+tag+'.pdf'
    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(1, 1)
    ax0 = plt.subplot(gs[0])
    ax0.plot(wlhr, Ahr, alpha=0.5, c='k')
    if telescope.mode != 'Imaging':
        ax0.plot(lam, A, alpha=0.7, color='orange', drawstyle='steps-mid', lw=2.0)
    else:
        ax0.plot(lam, A, 'o', alpha=0.7, color='orange', ms = 10.0)
        telescope.filter_wheel.plot(ax=ax0)
    ax0.errorbar(lam, spec, yerr=sig, fmt='o', color='k', ms=10.0)
    ax0.set_ylabel('Fp/Fs')
    ax0.set_xlabel('Wavelength [$\mu$m]')
    ax0.set_xlim([lammin,lammax])
    ax0.set_ylim([Amin, Amax])
    #ax0.set_ylim([-0.01,ymax+0.1])
    #ax0.set_ylim([-0.01,1.01])
    ax0.text(0.99, 0.99, plot_text,\
         verticalalignment='top', horizontalalignment='right',\
         transform=ax0.transAxes,\
         color='black', fontsize=20)
    # Save plot if saveplot==True
    if saveplot:
        fig.savefig(plot_tag)
        print 'Saved: '+plot_tag
    fig.show()

def plot_coronagraph_spectrum(wl, ofrat, sig, itime, d, ref_lam, SNR,
                              truth=None,
                              xlim=None, ylim=None,
                              title="",
                              save=False, tag=""):

    # Set matplotlib params
    mpl.rc('font', family='Times New Roman')
    mpl.rcParams['font.size'] = 25.0

    # Create figure
    fig = plt.figure(figsize=(12,10))
    gs = gridspec.GridSpec(1,1)
    ax = plt.subplot(gs[0])

    # Set string for plot text
    if itime > 2.0:
        timestr = "{:.0f}".format(itime)+' hours'
    else:
        timestr = "{:.0f}".format(itime*60)+' mins'
    plot_text = r'Distance = '+"{:.1f}".format(d)+' pc'+\
    '\n Integration time = '+timestr

    # If a reference wavelength is specified then return the SNR at that wl
    # corresponding to the integration time given
    if ref_lam:
        ireflam = find_nearest(wl,ref_lam)
        ref_SNR = SNR[ireflam]
        plot_text = plot_text + '\n SNR = '+"{:.1f}".format(ref_SNR)+\
            ' at '+"{:.2f}".format(wl[ireflam])+r' $\mu$m'

    if truth is not None:
        ax.plot(wl, truth*1e9, lw=2.0, color="purple", alpha=0.7, ls="steps-mid")
    ax.errorbar(wl, ofrat*1e9, yerr=sig*1e9, fmt='o', color='k', ms=5.0)

    ax.set_ylabel(r"F$_p$/F$_s$ ($\times 10^9$)")
    ax.set_xlabel("Wavelength [$\mu$m]")
    ax.set_title(title)
    ax.text(0.99, 0.99, plot_text,\
         verticalalignment='top', horizontalalignment='right',\
         transform=ax.transAxes,\
         color='black', fontsize=20)

    if ylim is not None: ax.set_ylim(ylim)
    if xlim is not None: ax.set_xlim(xlim)

    if save:
        fig.savefig(title+tag+".pdf")
    else:
        plt.show()

def process_noise(Dt, Cratio, cp, cb):
    """
    Computes SNR, noised data, and error on noised data.

    Parameters
    ----------
    Dt : float
        Telescope integration time in seconds
    Cratio : array
        Planet/Star flux ratio in each spectral bin
    cp : array
        Planet Photon count rate in each spectral bin
    cb : array
        Background Photon count rate in each spectral bin

    Returns
    -------
    cont : array
        Noised Planet/Star flux ratio in each spectral bin
    sigma : array
        One-sigma errors on flux ratio in each spectral bin
    SNR : array
        Signal-to-noise ratio in each spectral bin
    """

    # Calculate signal-to-noise assuming background subtraction (the "2")
    SNR  = cp*Dt/np.sqrt((cp + 2*cb)*Dt)

    # Calculate 1-sigma errors
    sigma= Cratio/SNR

    # Add gaussian noise to flux ratio
    cont = Cratio + np.random.randn(len(Cratio))*sigma

    return cont, sigma, SNR

def calc_SNR(itime, cp, cb, poisson=2.):

    cnoise = cp + poisson*cb
    SNR = cp*itime / np.sqrt(cnoise*itime)

    return SNR

def draw_noisy_spec(spectrum, SNR, apparent=False):

    if apparent:
        # Scale geometric albedo to apparent albedo (as if lambertian sphere) for quadrature
        scale = 3.*np.pi/2.
    else:
        scale = 1.0

    # Compute "sigma"
    sigma = scale * spectrum / SNR
    # Draw data points from normal distribution
    spec_noise = np.random.randn(len(spectrum))*sigma + scale*spectrum

    return spec_noise, sigma

def interp_cont_over_band(lam, cp, icont, iband):
    """
    Interpolate the continuum of a spectrum over a masked absorption or emission
    band.

    Parameters
    ----------
    lam : array
        Wavelength grid (abscissa)
    cp : array
        Planet photon count rates or any spectrum
    icont : list
        Indicies of continuum (neighboring points)
    iband : list
        Indicies of spectral feature (the band)

    Returns
    -------
    ccont : list
        Continuum planet photon count rates across spectral feature, where
        len(ccont) == len(iband)
    """
    # Linearly interpolate continuum points to band
    ccont = np.interp(lam[sorted(iband)], lam[sorted(icont)], cp[sorted(icont)])
    return ccont

def exptime_band(cp, ccont, cb, iband, SNR=5.0):
    """
    Calc the exposure time necessary to get a given S/N on a molecular band
    following Eqn 7 from Robinson et al. 2016.


    Parameters
    ----------
    cp :
        Planet count rate
    ccont :
        Continuum count rate
    cb :
        Background count rate
    iband :
        Indicies of molecular band
    SNR :
        Desired signal-to-noise ratio on molecular band

    Returns
    -------
    Telescope exposure time [hrs]
    """

    numerator = np.sum(cp[iband] + 2.*cb[iband])
    denominator = np.power(np.sum(np.fabs(ccont - cp[iband])),2)

    return np.power(SNR, 2) * numerator / denominator / 3600.0

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
