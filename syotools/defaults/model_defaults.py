#!/usr/bin/env python
"""
Created on Sat Oct 15 10:59:16 2016

@author: gkanarek, tumlinson
"""
from __future__ import (print_function, division, absolute_import, with_statement,
                        nested_scopes, generators)

import os 

#pathlib not supported in python 2
try:
    from pathlib import Path
    use_pathlib = True
except ImportError:
    import os
    use_pathlib = False

use_pathlib = False

from collections import OrderedDict
from syotools.utils import ordered_load

cwd = os.getenv('LUVOIR_SIMTOOLS_DIR')

#Load data from ascii table file (need a better method? maybe a FITS table?)
#and establish the default file path
if use_pathlib:
    spec_default_path = str(cwd+'syotools/'+'data/LUMOS_ETC.fits')
    pol_default_path = str(cwd+'syotools/'+'data/POLLUX_ETC.fits')
    yaml_default_path = str(cwd+'syotools/'+'defaults/model_defaults.yaml')
else:
    current_dir = os.path.dirname(os.path.abspath(__file__))
    spec_default_path = os.path.join(current_dir, '..', 'data', 'LUMOS_ETC.fits')
    pol_default_path = os.path.join(current_dir, '..', 'data', 'POLLUX_ETC.fits')
    yaml_default_path = os.path.join(current_dir, 'model_defaults.yaml')

#Now we load the defaults from model_defaults.yaml
#We use LUVOIR prelim values as defaults for the telescope & camera, and
#LUMOS prelim values as defaults for the spectrograph.
#Default exposure parameters are taken from the default HDI_ETC tool values
#   --> the _sed, _snr, and _magnitude default values are placeholders
#Default spectrograph parameters are taken from the default LUMOS_ETC tool values
#   --> the description, bef, R, wrange, wave, and aeff default values are placeholders
with open(yaml_default_path, 'r') as stream:
    all_defaults = ordered_load(stream)
    
default_telescope = all_defaults['Telescope']
default_camera = all_defaults['Camera']
default_exposure = all_defaults['Exposure']
default_spectrograph = OrderedDict([("_lumos_default_file", spec_default_path)])
default_spectrograph.update(all_defaults['Spectrograph'])
default_spectropolarimeter = OrderedDict([("_lumos_default_file", pol_default_path)])
default_spectropolarimeter.update(all_defaults['Spectropolarimeter'])
default_coronagraph = all_defaults['Coronagraph'] #placeholder

"""
##OLD VERSION FOR POSTERITY:

default_telescope = OrderedDict([(key, pre_encode(value)) for key,value in [
                                        ('name', 'LUVOIR'),
                                        ('aperture', 10. * u.m),
                                        ('temperature', 270. * u.K),
                                        ('ota_emissivity', 0.09 * u.dimensionless_unscaled),
                                        ('diff_limit_wavelength', 500. * u.nm)]])


#Again, prelim LUVOIR camera properties
default_camera = OrderedDict([(key, pre_encode(value)) for key,value in [
                                        ('name', 'HDI'),
                                        ('pivotwave', np.array([155., 228., 360., 440., 550., 640., 
                                                                790., 1260., 1600., 2220.]) * u.nm),
                                        ('bandnames', ['FUV', 'NUV', 'U','B','V','R','I', 
                                                       'J', 'H', 'K']),
                                        ('channels', [([0,1], 2), ([2, 3, 4, 5, 6], 2), ([7, 8, 9], 7)]),
                                        ('ab_zeropoint', np.array([35548., 24166., 15305., 12523., 
                                                                   10018., 8609., 6975., 4373., 
                                                                   3444., 2482.]) * (u.photon / u.s / u.cm**2 / u.nm)),
                                        ('total_qe', np.array([0.1, 0.1, 0.15, 0.45, 0.6, 0.6, 0.6, 
                                                               0.6, 0.6, 0.6]) * u.electron / u.ph),
                                        ('ap_corr', np.full(10, 1., dtype=float) * u.dimensionless_unscaled),
                                        ('bandpass_r', np.full(10, 5., dtype=float) * u.dimensionless_unscaled),
                                        ('dark_current', np.array([0.0005, 0.0005, 0.001, 0.001, 
                                                                   0.001, 0.001, 0.001, 0.002, 0.002, 
                                                                   0.002]) * (u.electron / u.s / u.pixel)),
                                        ('detector_rn', np.array([3., 3., 3., 3., 3., 3., 3., 4., 
                                                                  4., 4.]) * (u.electron / u.pixel))]])

#Default exposure parameters taken from HDI_ETC tool
default_exposure = {key: pre_encode(value) for key,value in [
                            ('_sed', np.full(10, 0.) * u.ABmag), #Set via sed_id by default, so this is a placeholder
                            ('_sed_id', 'fab'),
                            ('n_exp', 1),
                            ('_exptime', np.full(10, 1.) * u.hr),
                            ('_snr', np.full(10, 0.) * u.dimensionless_unscaled), #default unknown, so this is a placeholder
                            ('_magnitude', np.full(10, 0.) * u.ABmag), #Intepolated from SED by default, so this is a placeholder
                            ('unknown', 'snr')
                        ]}

#LUVOIR Multi-Object Spectrograph
 
default_spectrograph = {'name': 'LUMOS',
                        'modes': {'G120M': 'Med_Res_BEF', 
                                  'G150M': 'Med_Res_BEF', 
                                  'G180M': 'Med_Res_BEF', 
                                  'G155L': 'Low_Res_BEF', 
                                  'G145LL': 'LL_Mode_BEF'},
                        'befs': {mode: pre_encode(spec_default[mode+'_BEF'] * (u.erg / u.s / u.cm**3 / u.sr)) 
                                       for mode in ['Med_Res', 'Low_Res', 'LL_mode']},
                        'Rs': {'G120M': 30000., 
                               'G150M': 30000., 
                               'G180M': 30000., 
                               'G155L': 5000., 
                               'G145LL': 500.},
                        'ranges': {'G120M': pre_encode(np.array([1000., 1425.]) * u.AA), 
                                   'G150M': pre_encode(np.array([1225., 1600.]) * u.AA),  
                                   'G180M': pre_encode(np.array([1550., 1900.]) * u.AA), 
                                   'G155L': pre_encode(np.array([1000., 2000.]) * u.AA), 
                                   'G145LL': pre_encode(np.array([900., 1425.]) * u.AA)},
                        'mode': 'G150M',
                        'wave': pre_encode(spec_default['Wave'] * u.A),
                        'aeff': pre_encode(spec_default['A_eff'] * u.m**2)}

#Placeholder for default coronagraph

default_coronagraph = {}"""
