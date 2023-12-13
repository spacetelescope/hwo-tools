# Import dependent modules
import numpy as np
import sys
from .count_rates import count_rates
from .Noise import Output
import pdb

__all__ = ['count_rates_wrapper']

def count_rates_wrapper(Ahr, lamhr, solhr,
                   telescope, planet, star,
                   wantsnr=10.0, FIX_OWA = False, COMPUTE_LAM = False,
                   SILENT = False, NIR = True, THERMAL = False):
    """

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

    # Planet Parameters
    alpha  = planet.alpha           # Phase angle
    Phi    = planet.Phi             # Phase function
    Rp     = planet.Rp              # Planet radius (R_earth)
    r      = planet.a               # Semi-major axis (AU)
    d      = planet.distance        # Planet distance (pc)
    Nez    = planet.Nez             # Number of exo-zodi

    # Stellar Parameters
    Teff   = star.Teff              # Stellar effective temperature (K)
    Rs     = star.Rs                # Stellar radius (R_sun)

    # Telescope Parameters
    mode         = telescope.mode
    filter_wheel = telescope.filter_wheel
    lammin = telescope.lammin       # Wavelength minimum
    lammax = telescope.lammax       # Wavelength maximum
    Res    = telescope.resolution   # Resolving power
    diam   = telescope.diameter     # Diameter (m)
    Tput   = telescope.throughput   # Throughput
    C      = telescope.contrast     # Raw contrast
    IWA    = telescope.IWA          # Inner working angle
    OWA    = telescope.OWA          # Outer working angle
    Tsys   = telescope.Tsys         # Telescope/System temperature (K)
    Tdet   = telescope.Tdet         # Detector temperature (K)
    emis   = telescope.emissivity   # Emissivity

    # Set key system parameters
    De     = telescope.darkcurrent  # dark current (s**-1)
    DNHpix = telescope.DNHpix       # horizontal pixel spread of IFS spectrum
    Re     = telescope.readnoise    # read noise per pixel
    Dtmax  = telescope.Dtmax        # maximum exposure time (hr)
    X      = telescope.X            # size of photometric aperture (lambda/D)
    qe     = telescope.qe           # quantum efficiency

    # Set astrophysical parameters
    MzV  = 23.0                     # zodiacal light surface brightness (mag/arcsec**2)
    MezV = 22.0                     # exozodiacal light surface brightness (mag/arcsec**2)


    # Calculate count rates
    lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, DtSNR = \
        count_rates(Ahr, lamhr, solhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez,
                    mode   = mode,
                    filter_wheel = filter_wheel,
                    lammin = lammin,
                    lammax = lammax,
                    Res    = Res,
                    diam   = diam,
                    Tput   = Tput,
                    C      = C,
                    IWA    = IWA,
                    OWA    = OWA,
                    Tsys   = Tsys,
                    Tdet   = Tdet,
                    emis   = emis,
                    De     = De,
                    DNHpix = DNHpix,
                    Re     = Re,
                    Dtmax  = Dtmax,
                    X      = X,
                    qe     = qe,
                    MzV    = MzV,
                    MezV   = MezV,
                    wantsnr=wantsnr, FIX_OWA=FIX_OWA, COMPUTE_LAM=COMPUTE_LAM,
                    SILENT=SILENT, NIR=NIR, THERMAL=THERMAL)

    # Cram all the coronagraph output arrays into an Output object
    output = Output(lam=lam, dlam=dlam, A=A, q=q, Cratio=Cratio,
                    cp=cp, csp=csp, cz=cz, cez=cez, cD=cD, cR=cR,
                    cth=cth, DtSNR=DtSNR)

    # Return output object
    return output
