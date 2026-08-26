
#!/usr/bin/env python
# coding: utf-8

# ## Basic Operations of the SYOTools version of the HWO ETCs

#first the necessary imports 
from syotools.models import Camera, Telescope, Spectrograph, SourcePhotometricExposure, SourceSpectrographicExposure
from syotools.utils.jsonunit import str_jsunit

# create a Telescope and Exposure 
# and select a Camera from the Telescope
e, t = SourcePhotometricExposure(), Telescope()
t.set_from_hwome('EAC5')
c = t.instruments["HRI_S.HRI_S_UVIS_Imager"]

# if you do no further customization to the exposure, you get the defaults
c.add_exposure(e)

for attr in t._tracked_attributes:
    print('{}: {}'.format(attr, str_jsunit(getattr(t,attr))))

# set a camera band
c.band = "HRI_S_UVIS.HRI_Johnson_V"

# set what you're calculating for (this also kicks off that kind of calculation)
e.unknown="snr"

# print the answer
print(e.snr)

