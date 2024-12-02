#!/usr/bin/env python
# coding: utf-8

# ## Basic Operations of the SYOTools version of the HWO ETCs

from syotools.models import Camera, Telescope, Spectrograph, Source, SourcePhotometricExposure
from syotools.utils.jsonunit import str_jsunit
import numpy as np, astropy.units as u 

# create a Telescope, Camera, and Exposure 
# if called without arguments, you get the defaults 
t, c = Telescope(), Camera()

#set the telescope to EAC1 
t.set_from_json('EAC1')

source = Source() 
source.set_sed('Flat (AB)', 30., 0.0, 0.0, 'galex,fuv')   

e = SourcePhotometricExposure() 
e.source = source

t.add_camera(c)
c.add_exposure(e)

for attr in t._tracked_attributes:
    print('{}: {}'.format(attr, str_jsunit(getattr(t,attr))))


for attr in c._tracked_attributes:
    print('{}: {}'.format(attr, str_jsunit(getattr(c,attr))))


for attr in e._tracked_attributes:
    print('{}: {}'.format(attr, str_jsunit(getattr(e, attr))))

print() 
print() 
print() 
print('----------- Setting SNR as the Unknown -----------') 
print(e.snr)

print() 
print() 
print() 
print('----------- Setting Exptime as the Unknown -----------') 
source.set_sed('Flat (AB)', 30., 0.0, 0.0, 'galex,fuv')   
e.source = source
e._snr = np.array([10, 10, 10, 10, 10, 10, 10, 10, 10, 10]) * u.Unit('electron(1/2)') 
e.unknown = 'exptime' 
print(e.exptime) 

print() 
print() 
print() 
print('----------- Setting Magnitude as the Unknown -----------') 
source.set_sed('Flat (AB)', 30., 0.0, 0.0, 'galex,fuv')   
e.source = source
e._snr = np.array([10, 10, 10, 10, 10, 10, 10, 10, 10, 10]) * u.Unit('electron(1/2)') 
#e._exptime = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) * u.Unit('hr') 
e._exptime = [np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1]), 'hr'] 
e.unknown = 'magnitude' 
print(e.magnitude) 




