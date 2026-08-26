#!/usr/bin/env python
# coding: utf-8

# ## Basic Operations of the SYOTools version of the HWO ETCs

from syotools.models import Camera, Telescope, Spectrograph, Source, SourcePhotometricExposure
from syotools.utils.jsonunit import str_jsunit
import numpy as np, astropy.units as u 

# create a Telescope
t = Telescope()

#set the telescope to EAC5 with HRI_S.HRI_S_UVIS_Imager
t.set_from_hwome("EAC5")
c = t.instruments["HRI_S.HRI_S_UVIS_Imager"]

# Create a Source
source = Source() 
source.set_sed('Flat (AB)', 30., 0.0, 0.0, 'galex,fuv')   

# Create an Exposure, and add the Source to the Exposure
e = SourcePhotometricExposure() 
e.source = source

# Add the Exposure to the Camera
c.add_exposure(e)

# Set the exposure time of the Exposure
e.exptime = 1 * u.hr # anything Astropy can convert to seconds is fine

for attr in t._tracked_attributes:
    print('{}: {}'.format(attr, str_jsunit(getattr(t,attr))))


print() 
print() 
print() 
print('----------- Setting SNR as the Unknown -----------') 
# Setting the unknown needs to be the final step of setting up the exposure.
# It unlocks auto-recalc - after this point, SYOTools will recompute any time
# input data changes so that the values in its arrays are correct when queried.
e.unknown = "snr" 
print(e.snr)

print() 
print() 
print() 
print('----------- Setting Exptime as the Unknown -----------') 
# set this first, otherwise attempting to set the SNR will fail (because without this it's still listed as the unknown value to solve for)
e.unknown = 'exptime'
e.snr = 10
# notice we don't reset the source. That carries over.
print(e.exptime) 

print() 
print() 
print() 
print('----------- Setting Magnitude as the Unknown -----------') 
source.set_sed('Flat (AB)', 30., 0.0, 0.0, 'galex,fuv')   
e.source = source
e.snr = 10 * u.dimensionless_unscaled # also a valid way to consider SNR
e.exptime = 1 * u.hr
e.unknown = 'magnitude'
print(e.magnitude)
# Notice that nothing happens when we set a new source.
source.set_sed('G2V Star', 28., 0.0, 0.0, 'galex,fuv')
print(e.magnitude)
# We can make the system run a new calculation. That will actually change numbers.
e.calculate_magnitude()
print(e.magnitude)
# Where the Source requires manual intervention to recompute, we can reset the snr
# and SYOTools will recalculate automatically.
e.snr = 50 * u.dimensionless_unscaled
print(e.magnitude)




