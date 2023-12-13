"""
Telescope, Planet, and Star classes.
"""

import numpy as np
from .utils import Input, Loadin

__all__ = ['Telescope', 'Planet', 'Star']

################################################################################
# TELESCOPE
################################################################################

class Telescope(object):
    """
    A class to represent a telescope object and all design specifications
    therein

    Parameters
    ----------
    mode : str
        Telescope observing modes: 'IFS', 'Imaging'
    lammin : float
        Minimum wavelength (um)
    lammax : float
        Maximum wavelength (um)
    R : float
        Spectral resolution (lambda / delta-lambda)
    Tsys : float
        Telescope temperature (K)
    D : float
        Telescope diameter (m)
    emis : float
        Telescope emissivity
    IWA : float
        Inner Working Angle (lambda/D)
    OWA : float
        Outer Working Angle (lambda/D)
    Tput : float
        Telescope throughput
    C : float
        Raw Contrast
    De : float
        Dark current (s**-1)
    DNHpix : float
        Horizontal pixel spread of IFS spectrum
    Re : float
        Read noise per pixel
    Dtmax : float
        Maximum exposure time (hr)
    X : float
        Size of photometric aperture (lambda/D)
    q : float
        Quantum efficiency
    filter_wheel : Wheel (optional)
        Wheel object containing imaging filters

    Methods
    -------
    from_file(path)
        Initialize telescope object using telescope parameters in the Input file
    default_luvoir()
        Initialize telescope object using current LUVOIR parameters (Not decided!)
    default_habex()
        Initialize telescope object using current HabEx parameters (Not decided!)
    default_wfirst()
        Initialize telescope object using current WFIRST parameters (Not decided!)
    """

    # Define a constructor
    def __init__(self, mode='IFS', lammin=0.3,lammax=2.0,R=70.,Tput=0.2,\
                 D=8.0,Tsys=274.,Tdet=50.,IWA=0.5, OWA=30000.,emis=0.9,\
                 C=1e-10,De=1e-4,DNHpix=3.,Re=0.1,Dtmax=1.0,X=0.7,q=0.9,\
                 filter_wheel=None):
        self._mode=mode
        self.lammin=lammin
        self.lammax=lammax
        self.resolution=R
        self.throughput=Tput
        self.diameter=D
        self.Tsys=Tsys
        self.Tdet=Tdet
        self.IWA=IWA
        self.OWA=OWA
        self.emissivity=emis
        self.contrast=C

        self.darkcurrent=De
        self.DNHpix=DNHpix
        self.readnoise=Re
        self.Dtmax=Dtmax
        self.X=X
        self.qe=q

        self._filter_wheel=filter_wheel

        if self._mode == 'Imaging':
            from filters.imager import johnson_cousins
            self._filter_wheel = johnson_cousins()

    @classmethod
    def from_file(cls, path):

        # Read-in Telescope params using Loadin class
        L = Loadin(path)

        # Return new class instance
        return cls(mode=L.mode, lammin=L.lammin, lammax=L.lammax, R=L.resolution,
                   Tput=L.throughput, D=L.diameter, Tsys=L.Tsys, Tdet=L.Tdet, IWA=L.IWA,
                   OWA=L.OWA, emis=L.emissivity, C=L.contrast, De=L.darkcurrent,
                   DNHpix=L.DNHpix, Re=L.readnoise, Dtmax=L.Dtmax, X=L.X,
                   q=L.qe, filter_wheel=L.filter_wheel)

    @classmethod
    def default_luvoir(cls):
        # Return new class instance
        return cls(mode="IFS", lammin=0.5, lammax=1.0, R=70.,
                   Tput=0.05, D=12., Tsys=150., Tdet=50., IWA=3.0,
                   OWA=20.0, emis=0.9, C=1e-10, De=1e-4,
                   DNHpix=3.0, Re=0.1, Dtmax=1.0, X=1.5,
                   q=0.9, filter_wheel=None)

    @classmethod
    def default_habex(cls):
        print "These HabEx parameters are not confirmed yet!"
        # Return new class instance
        return cls(mode="IFS", lammin=0.4, lammax=2.5, R=70.,
                   Tput=0.05, D=6., Tsys=150., Tdet=50., IWA=3.0,
                   OWA=20.0, emis=0.9, C=1e-10, De=1e-4,
                   DNHpix=3.0, Re=0.1, Dtmax=1.0, X=1.5,
                   q=0.9, filter_wheel=None)

    @classmethod
    def default_wfirst(cls):
        print "These WFIRST parameters are not confirmed yet!"
        # Return new class instance
        return cls(mode="IFS", lammin=0.6, lammax=1.0, R=70.,
                   Tput=0.05, D=2.4, Tsys=150., Tdet=50., IWA=3.0,
                   OWA=20.0, emis=0.9, C=1e-9, De=1e-4,
                   DNHpix=3.0, Re=0.1, Dtmax=1.0, X=1.5,
                   q=0.9, filter_wheel=None)

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, value):
        self._mode = value
        if value == 'Imaging':
            from filters.imager import johnson_cousins
            self._filter_wheel = johnson_cousins()
        else:
            self._filter_wheel = None

    @property
    def filter_wheel(self):
        return self._filter_wheel

    @filter_wheel.setter
    def filter_wheel(self, value):
        if (value.__class__.__name__ == 'Wheel') or (value.__class__.__base__.__name__ == 'Wheel'):
            self._filter_wheel = value
        else:
            print "Error in Telescope: Specified filter wheel is not of type 'Wheel'"
            self._filter_wheel = None

    def __str__(self):
        string = 'Coronagraph: \n------------\n'+\
            '- Telescope observing mode : '+"%s" % (self.mode)+'\n'+\
            '- Minimum wavelength (um) : '+"%s" % (self.lammin)+'\n'+\
            '- Maximum wavelength (um)  : '+"%s" % (self.lammax)+'\n'+\
            '- Spectral resolution (lambda / delta-lambda)  : '+"%s" % (self.resolution)+' \n'+\
            '- Telescope/System temperature (K)  : '+"%s" % (self.Tsys)+' \n'+\
            '- Detector temperature (K)  : '+"%s" % (self.Tdet)+' \n'+\
            '- Telescope diameter (m)  : '+"%s" % (self.diameter)+' \n'+\
            '- Telescope emissivity  : '+"%s" % (self.emissivity)+' \n'+\
            '- Inner Working Angle (lambda/D)  : '+"%s" % (self.IWA)+' \n'+\
            '- Outer Working Angle (lambda/D)  : '+"%s" % (self.OWA)+' \n'+\
            '- Telescope throughput  : '+"%s" % (self.throughput)+' \n'+\
            '- Raw Contrast  : '+"%s" % (self.contrast)+' \n'+\
            '- Dark current (s**-1)  : '+"%s" % (self.darkcurrent)+' \n'+\
            '- Horizontal pixel spread of IFS spectrum  : '+"%s" % (self.DNHpix)+' \n'+\
            '- Read noise per pixel  : '+"%s" % (self.readnoise)+' \n'+\
            '- Maximum exposure time (hr)  : '+"%s" % (self.Dtmax)+' \n'+\
            '- Size of photometric aperture (lambda/D)  : '+"%s" % (self.X)+' \n'+\
            '- Quantum efficiency  : '+"%s" % (self.qe)+''
        return string

################################################################################
# PLANET
################################################################################

def lambertPhaseFunction(alpha):
    """
    Calculate the Lambertian Phase Function from the phase angle.

    Parameters
    ----------
    alpha: float
        Planet phase angle (degrees)

    Returns
    -------
    Phi : float
        Lambertian phase function
    """
    alpha = alpha * np.pi / 180.
    print 'phi = ' +str((np.sin(alpha) + (np.pi - alpha) * np.cos(alpha)) / np.pi)
    return (np.sin(alpha) + (np.pi - alpha) * np.cos(alpha)) / np.pi

class Planet(object):
    """
    A class to represent a planet and all associated parameters of the planet
    to be observed.

    Parameters
    ----------
    name : string
        Planet name from database
    star : string
        Stellar type of planet host star
    d : float
        Distance to system (pc)
    Nez : float
        Number of exzodis (zodis)
    Rp : float
        Radius of planet (Earth Radii)
    a : float
        Semi-major axis (AU)
    alpha : float
        Phase angle (deg)
    Phi : float
        Lambertian phase function
    MzV : float
        Zodiacal light surface brightness (mag/arcsec**2)
    MezV : float
        exozodiacal light surface brightness (mag/arcsec**2)

    Methods
    -------
    from_file()
        Initialize object using planet parameters in the Input file
    """

    # Define a constructor
    def __init__(self, name='earth', star='sun', d=10.0,Nez=1.0,\
                 Rp=1.0, a=1.0, alpha=90.,\
                 MzV=23.0, MezV=22.0):
        self.name=name
        self.star=star
        self.distance=d
        self.Nez=Nez
        self.Rp=Rp
        self.a=a
        self._alpha=alpha
        self._Phi = None
        self.MzV  = MzV     # zodiacal light surface brightness (mag/arcsec**2)
        self.MezV = MezV     # exozodiacal light surface brightness (mag/arcsec**2)

        if self._Phi is None:
            self._Phi = lambertPhaseFunction(self._alpha)
        else:
            raise Exception("Error in Planet Phase Function (Phi)")

    @classmethod
    def from_file(cls, path):

        # Read-in Telescope params using Loadin class
        L = Loadin(path)

        # Return new class instance
        return cls(name=L.name, star=L.star, d=L.distance, Nez=L.Nez,
                   Rp=L.Rp, a=L.a, alpha=L.alpha, MzV=L.MzV, MezV=L.MezV)

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = value
        self._Phi = lambertPhaseFunction(value)

    @property
    def Phi(self):
        return self._Phi

    @Phi.setter
    def Phi(self, value):
        self._Phi = value

    def __str__(self):
        string = 'Planet: \n-------\n'+\
            '- Planet name  : '+"%s" % (self.name)+'\n'+\
            '- Stellar type of planet host star : '+"%s" % (self.star)+'\n'+\
            '- Distance to system (pc) : '+"%s" % (self.distance)+'\n'+\
            '- Number of exzodis (zodis) : '+"%s" % (self.Nez)+'\n'+\
            '- Radius of planet (Earth Radii) : '+"%s" % (self.Rp)+'\n'+\
            '- Semi-major axis (AU) : '+"%s" % (self.a)+'\n'+\
            '- Phase angle (deg) : '+"%s" % (self.alpha)+'\n'+\
            '- Lambertian phase function : '+"%s" % (self.Phi)+'\n'+\
            '- Zodiacal light surface brightness (mag/arcsec**2) : '+"%s" % (self.MzV)+'\n'+\
            '- Exozodiacal light surface brightness (mag/arcsec**2) : '+"%s" % (self.MezV)
        return string

################################################################################
# STAR
################################################################################

class Star(object):
    """
    A class to represent the stellar host for an exoplanet observation

    Parameters
    ----------
    Teff : float
        Stellar effective temperature [K]
    Rs : float
        Stellar radius [Solar Radii]
    """

    def __init__(self, Teff=5780.0, Rs=1.0):
        self.Teff=Teff
        self.Rs=Rs

    @classmethod
    def from_file(cls, path):

        # Read-in Telescope params using Loadin class
        L = Loadin(path)

        # Return new class instance
        return cls(Teff=L.Teff, Rs=L.Rs)

    def __str__(self):
        string = 'Star: \n-----\n'+\
            '- Effective Temperature (K) : '+"%s" % (self.Teff)+'\n'+\
            '- Radius (Solar Radii) : '+"%s" % (self.Rs)
        return string
