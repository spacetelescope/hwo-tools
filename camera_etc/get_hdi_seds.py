from __future__ import print_function
import pysynphot as S 
import os 

def add_spectrum_to_library():

    # cwd = os.getenv('LUVOIR_SIMTOOLS_DIR')

    spec_dict = {}

    flatsp = S.FlatSpectrum(30, fluxunits='abmag') 
    flatsp = flatsp.renorm(30., 'abmag', S.ObsBandpass('johnson,v'))
    flatsp.convert('abmag') 
    flatsp.convert('nm') 
    spec_dict['Flat (AB)'] = flatsp 

    bb = S.BlackBody(5000)
    bb.convert('abmag') 
    bb.convert('nm') 
    spec_dict['Blackbody (5000K)'] = bb 

    filename_qso = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'qso_template.fits')
    qso = S.FileSpectrum(filename_qso)
    q = qso.renorm(30., 'abmag', S.ObsBandpass('galex,fuv'))
    q.convert('nm') 
    q.convert('abmag') 
    spec_dict['QSO'] = q

    filename_star = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'pickles', 'dat_uvk', 'pickles_uk_1.fits')
    star = S.FileSpectrum(filename_star)
    s = star.renorm(30., 'abmag', S.ObsBandpass('johnson,v'))
    s.convert('abmag') 
    s.convert('nm') 
    spec_dict['O5V Star'] = s

    filename_star = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'pickles', 'dat_uvk', 'pickles_uk_6.fits')
    star = S.FileSpectrum(filename_star)
    s = star.renorm(30., 'abmag', S.ObsBandpass('johnson,v'))
    s.convert('abmag') 
    s.convert('nm') 
    spec_dict['B5V Star'] = s

    filename_star = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'pickles', 'dat_uvk', 'pickles_uk_26.fits')
    star = S.FileSpectrum(filename_star)
    s = star.renorm(30., 'abmag', S.ObsBandpass('johnson,v'))
    s.convert('abmag') 
    s.convert('nm') 
    spec_dict['G2V Star'] = s

    filename_star = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'pickles', 'dat_uvk', 'pickles_uk_40.fits')
    star = S.FileSpectrum(filename_star)
    s = star.renorm(30., 'abmag', S.ObsBandpass('johnson,v'))
    s.convert('abmag') 
    s.convert('nm') 
    spec_dict['M2V Star'] = s

    filename = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'ngc1068_template.fits')
    agn  = S.FileSpectrum(filename)
    a = agn.renorm(30., 'abmag', S.ObsBandpass('johnson,v'))
    a.convert('abmag') 
    a.convert('nm') 
    spec_dict['NGC 1068'] = a

    filename = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'etc_models', 'el_cww_fuv_001.fits')
    ell  = S.FileSpectrum(filename)
    e = ell.renorm(30., 'abmag', S.ObsBandpass('johnson,v'))
    e.convert('abmag') 
    e.convert('nm') 
    spec_dict['Elliptical Galaxy'] = e

    filename = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'etc_models', 'sbc_cb2004a_001.fits')
    ell  = S.FileSpectrum(filename)
    e = ell.renorm(30., 'abmag', S.ObsBandpass('johnson,v'))
    e.convert('abmag') 
    e.convert('nm') 
    spec_dict['Sbc Galaxy'] = e

    filename = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'etc_models', 'sb1_kinney_fuv_001.fits')
    ell  = S.FileSpectrum(filename)
    e = ell.renorm(30., 'abmag', S.ObsBandpass('johnson,v'))
    e.convert('abmag') 
    e.convert('nm') 
    spec_dict['Starburst Galaxy'] = e

    filename = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'galactic', 'orion_template.fits')
    neb  = S.FileSpectrum(filename)
    n = neb.renorm(30., 'abmag', S.ObsBandpass('johnson,v'))
    n.convert('abmag') 
    n.convert('nm') 
    spec_dict['Orion Nebula'] = n

    return spec_dict

