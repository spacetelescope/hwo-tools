#!/usr/bin/env python
"""
Created on Tue Oct 18 14:10:45 2016

@author: gkanarek
"""

from __future__ import (print_function, division, absolute_import, with_statement,
                        nested_scopes, generators)

import os

import syotools.cdbs

import pysynphot as pys
import astropy.io.ascii as asc
from syotools.utils import pre_encode

#pathlib is not supported in python 2
try:
    from pathlib import Path
    use_pathlib = True
except ImportError:
    use_pathlib = False

pysyn_base = os.environ['PYSYN_CDBS']
data_base = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data'))

def load_txtfile(spec):
    fname = spec['file']
    band = spec.get('band', 'johnson,v')
    if use_pathlib:
        path = Path(fname[0])
        for f in fname[1:]:
            path = path / f
        abspath = str(path.resolve())
    else:
        abspath = os.path.abspath(os.path.join(*fname))
    tab = asc.read(abspath, names=['wave','flux']) 
    sp = pys.ArraySpectrum(wave=tab['wave'], flux=tab['flux'], waveunits='Angstrom', fluxunits='flam')
    sp = sp.renorm(30., 'abmag', pys.ObsBandpass(band))
    sp.convert('abmag')
    sp.convert('nm')
    return pre_encode(sp)

def load_fesc(spec):
    fname = spec['file']
    band = spec.get('band', 'johnson,v')
    if use_pathlib:
        path = Path(fname[0])
        for f in fname[1:]:
            path = path / f
        abspath = str(path.resolve())
    else:
        abspath = os.path.abspath(os.path.join(*fname))
    tab = asc.read(abspath)
    sp = pys.ArraySpectrum(wave=tab['lam'], flux=tab['lh1=17.5'], 
                           waveunits='Angstrom', fluxunits='flam')
    sp = sp.renorm(30., 'abmag', pys.ObsBandpass(band))
    sp.convert('abmag')
    sp.convert('nm')

def load_pysfits(spec):
    fname = spec['file']
    band = spec.get('band', 'johnson,v')
    if use_pathlib:
        path = Path(fname[0])
        for f in fname[1:]:
            path = path / f
        abspath = str(path.resolve())
    else:
       abspath = os.path.abspath(os.path.join(*fname)) 
    sp = pys.FileSpectrum(abspath)
    sp = sp.renorm(30., 'abmag', pys.ObsBandpass(band))
    sp.convert('abmag')
    sp.convert('nm')
    return pre_encode(sp)

specs = {'ctts': {'desc': 'Classical T-Tauri Star', 
                  'file': [data_base, 'CTTS_etc_d140pc_101116.txt'],
                  'band': 'galex,fuv'},
         'ctts2': {'desc': 'Classical T-Tauri Star', 
                  'file': [data_base, 'ttauri.txt'],
                  'band': 'galex,fuv'},
         'mdwarf': {'desc': 'M1 Dwarf', 
                    'file': [data_base, 'dM1_etc_d5pc_101116.txt'],
                    'band': 'galex,fuv'},
         'mdwarf2': {'desc': 'M Dwarf', 
                    'file': [data_base, 'dM_d1pc_pollux.txt'],
                    'band': 'galex,fuv'},
         'gdwarf': {'desc': 'G Dwarf', 
                    'file': [data_base, 'dG_d5pc_pollux.txt'],
                    'band': 'galex,fuv'},
         's99': {'desc': '10 Myr Starburst', 
                 'file': [data_base, '10Myr_Starburst_nodust.dat'],
                 'band': 'galex,fuv'},
         'qso': {'desc': 'QSO', 
                 'file': [pysyn_base, 'grid', 'agn', 'qso_template.fits'], 
                 'band': 'galex,fuv'},
         'syfrt1': {'desc': 'Seyfert 1', 
                 'file': [data_base, 'Seyfert_1_template.txt'], 
                 'band': 'galex,fuv'},
         'syfrt2': {'desc': 'Seyfert 2', 
                 'file': [data_base,'Seyfert_2_template.txt'], 
                 'band': 'galex,fuv'},
         'liner': {'desc': 'Liner', 
                 'file': [data_base, 'LINER_template.txt'], 
                 'band': 'galex,fuv'},
         'o5v': {'desc': 'O5V Star', 
                 'file': [pysyn_base, 'grid', 'pickles', 'dat_uvk', 
                          'pickles_uk_1.fits']},
         'g2v': {'desc': 'G2V Star', 
                 'file': [pysyn_base, 'grid', 'pickles', 'dat_uvk', 
                          'pickles_uk_26.fits']},
         'orion': {'desc': 'Orion Nebula', 
                   'file': [pysyn_base, 'grid', 'galactic', 
                            'orion_template.fits']},
         'g191b2b': {'desc': 'G191B2B (WD)', 
                     'file': [pysyn_base, 'calspec', 'g191b2b_mod_010.fits']},
         'gd71': {'desc': 'GD71 (WD)',
                  'file': [pysyn_base, 'calspec', 'gd71_fos_003.fits'],
                  'band': 'galex,fuv'},
         'gd153': {'desc': 'GD153 (WD)',
                   'file': [pysyn_base, 'calspec', 'gd153_fos_003.fits'],
                   'band': 'galex,fuv'},
         'nodust': {'desc': 'Starburst, No Dust', 
                    'file': [pysyn_base, 'grid', 'kc96', 
                             'starb1_template.fits']},
         'ebv6': {'desc': 'Starburst, E(B-V) = 0.6', 
                  'file': [pysyn_base, 'grid', 'kc96', 
                           'starb6_template.fits']},
         'b5v': {'desc': 'B5V Star', 
                 'file': [pysyn_base, 'grid', 'pickles', 'dat_uvk', 
                          'pickles_uk_6.fits']},
         'm2v': {'desc': 'M2V Star', 
                 'file': [pysyn_base, 'grid', 'pickles', 'dat_uvk', 
                          'pickles_uk_40.fits']},
         'elliptical': {'desc': 'Elliptical Galaxy', 
                        #'file': [pysyn_base, 'grid', 'etc_models', 
                        #         'el_cww_fuv_001.fits']},
                        'file': [pysyn_base, 'grid', 'kc96', 
                                'elliptical_template.fits'],
                        'band': 'galex,fuv'},
         'sbc': {'desc': 'Sbc Galaxy', 
                 #'file': [pysyn_base, 'grid', 'etc_models', 
                 #         'sbc_cb2004a_001.fits']},
                 'file': [pysyn_base, 'grid', 'kc96', 
                          'sc_template.fits'],
                 'band': 'galex,fuv'},
         'starburst': {'desc': 'Starburst Galaxy',
                       #'file': [pysyn_base, 'grid', 'etc_models', 
                       #         'sb1_kinney_fuv_001.fits']},
                       'file': [pysyn_base, 'grid', 'kc96', 
                                'starb1_template.fits'],
                       'band': 'galex,fuv'},
         'ngc1068': {'desc': 'NGC 1068',
                     'file': [pysyn_base, 'grid', 'agn', 
                              'ngc1068_template.fits']},
         'hi1hei1': {'desc': 'Galaxy with f_esc, HI=1, HeI=1',
                     'file': [data_base, 'fesc', 'fe_lyccontrot1.000hi1.000hei.txt'],
                     'band': 'galex,fuv'},
         'hi0hei1': {'desc': 'Galaxy with f_esc, HI=0.001, HeI=1',
                     'file': [data_base, 'fesc', 'fe_lyccontrot0.001hi1.000hei.txt'],
                     'band': 'galex,fuv'}
        }

default_spectra = {'specs':{}, 'descs':{}}

for (specid, spec) in specs.items():
    default_spectra['descs'][specid] = spec['desc']
    fesc = 'fesc' in spec['file']
    fits = spec['file'][-1].endswith('fits')
    if fits:
        default_spectra['specs'][specid] = load_pysfits(spec)
    elif fesc:
        default_spectra['specs'][specid] = load_fesc(spec)
    else:
        default_spectra['specs'][specid] = load_txtfile(spec)

flatsp = pys.FlatSpectrum(30.0, fluxunits='abmag') 
flatsp = flatsp.renorm(30.0, 'abmag', pys.ObsBandpass('johnson,v'))
flatsp.convert('abmag') 
flatsp.convert('nm') 
default_spectra['specs']['fab'] = pre_encode(flatsp)
default_spectra['descs']['fab'] = 'Flat (AB)'

flamsp = pys.FlatSpectrum(30.0, fluxunits='flam')
flamsp = flamsp.renorm(30.0, 'abmag', pys.ObsBandpass('galex,fuv'))
flamsp.convert('abmag') 
flamsp.convert('nm') 
default_spectra['specs']['flam'] = pre_encode(flamsp)
default_spectra['descs']['flam'] = 'Flat in F_lambda'

bb = pys.BlackBody(5000)
bb.convert('abmag') 
bb.convert('nm') 
default_spectra['specs']['bb'] = pre_encode(bb)
default_spectra['descs']['bb'] = 'Blackbody (5000K)'
