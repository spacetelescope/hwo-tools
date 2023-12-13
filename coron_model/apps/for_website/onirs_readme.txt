The O/NIRS simulator tool is based off Tyler Robinson's coronagraph noise model (Robinson et al. 2016). The python version of this tool was produced by Jacob Lustig-Yaeger (available at https://github.com/jlustigy/coronagraph), and it was modified to simulate O/NIRS observations by Giada Arney by removing the coronagraphic component. The Bokeh rendering was based off of the coronagraph noise model tool developed by Jason Tumlinson and Giada Arney.

The values of internal (i.e. non-slider) input parameters used here are:

Dark current = 1e-4 counts/s
Read noise = 0.1 counts per pixel 
Detector temperature = 50 K
Mirror Emissivity = 0.09
DNHpix = 3.0 (= number of horizontal/spatial pixels in dispersed spectrum)
width of photometric aperture = 1.5 (times lambda/D)
quantum efficiency factor = 0.9 
MzV = 23.0 ( = V-band zodiacal light surface brightness)
Throughput = 15%

Note that this simulator differs from the simulator described in Robinson et al. (2016) because here we assume a constant quantum efficiency as a function of wavelength ( = 0.9), whereas Robinson et al. (2016) assumes wavelength dependence.

------------
