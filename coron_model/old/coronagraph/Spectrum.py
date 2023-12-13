import numpy as np
import readsmart as rs

class Spectrum(object):

    def __init__(self, wl=None, spec=None, starflux=None):
        self.wl = wl
        self.spec = spec
        self.starflux = starflux

    def plot(self):
        return

    @classmethod
    def from_smart(cls, path):
        # Check path extension
        if not path.endswith('.rad'):
            print("Incompatible input file.")
            return

        # Read-in .rad file with readsmart
        wlhr, wno, solar_spec, TOA_flux, rad_streams = readsmart.rad(radfile,getdata=True)

        # Calculate Hi-res reflectivity spectrum
        Ahr = (TOA_flux / solar_spec) #* np.pi / planet.Phi

        # Compose new class instance
        return cls(wl=wlhr, spec=Ahr, starflux=solar_spec)

    @classmethod
    def from_solar_system(cls, name=None, planetdir='planets/'):
        return
