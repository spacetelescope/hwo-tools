import pysynphot as S 
import os 
import numpy as np 
from astropy.io import ascii

def add_spectrum_to_library():

    cwd = os.getenv('LUVOIR_SIMTOOLS_DIR')
    cwd = './' 

    spec_dict = {}

    tab = ascii.read(cwd+'data/CTTS_etc_d140pc_101116.txt', names=['wave','flux'])
    sp = S.ArraySpectrum(wave=tab['wave'], flux=tab['flux'], waveunits='Angstrom', fluxunits='flam')
    ctts = sp.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['Classical T Tauri'] = ctts

    tab = ascii.read(cwd+'data/dM1_etc_d5pc_101116.txt', names=['wave','flux'])
    sp = S.ArraySpectrum(wave=tab['wave'], flux=tab['flux'], waveunits='Angstrom', fluxunits='flam')
    Mdwarf = sp.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['M1 Dwarf'] = Mdwarf

    tab = ascii.read(cwd+'data/10Myr_Starburst_nodust.dat', names=['wave', 'flux'])
    sp = S.ArraySpectrum(wave=tab['wave'], flux=tab['flux'], waveunits='Angstrom', fluxunits='flam')
    s99 = sp.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['10 Myr Starburst'] = s99

    filename_qso = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'qso_template.fits')
    qso = S.FileSpectrum(filename_qso)
    q = qso.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['QSO'] = q

    filename_star = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'pickles', 'dat_uvk', 'pickles_uk_1.fits')
    star = S.FileSpectrum(filename_star)
    s = star.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['O5V Star'] = s

    filename_star = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'pickles', 'dat_uvk', 'pickles_uk_26.fits')
    star = S.FileSpectrum(filename_star)
    s = star.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['G2V Star'] = s

    filename_star = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'galactic', 'orion_template.fits')
    star = S.FileSpectrum(filename_star)
    s = star.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['Orion Nebula'] = s

    filename_star = os.path.join(os.environ['PYSYN_CDBS'], 'calspec', 'g191b2b_mod_010.fits')
    star = S.FileSpectrum(filename_star)
    star_wave_sample = np.arange(2500) + 500. 
    star_flux_sample = star.sample(star_wave_sample) 
    sp = S.ArraySpectrum(wave=star_wave_sample, flux=star_flux_sample, waveunits='Angstrom', fluxunits='flam')
    s = sp.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['G191B2B (WD)'] = s

    filename_star = os.path.join(os.environ['PYSYN_CDBS'], 'calspec', 'gd71_fos_003.fits')
    star = S.FileSpectrum(filename_star)
    s = star.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['GD71 (WD)'] = s

    filename_star = os.path.join(os.environ['PYSYN_CDBS'], 'calspec', 'gd153_fos_003.fits')
    star = S.FileSpectrum(filename_star)
    s = star.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['GD153 (WD)'] = s

    filename_star = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'kc96', 'starb1_template.fits')
    star = S.FileSpectrum(filename_star)
    s = star.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['Starburst, No Dust'] = s

    filename_star = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'kc96', 'starb6_template.fits')
    star = S.FileSpectrum(filename_star)
    s = star.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['Starburst, E(B-V) = 0.6'] = s

    filename_fesc = cwd+'data/fesc/fe_lyccontrot1.000hi1.000hei.txt'
    tab = ascii.read(filename_fesc)
    sp = S.ArraySpectrum(wave=tab['lam'], flux=tab['lh1=17.5'], waveunits='Angstrom', fluxunits='flam')
    fesc1 = sp.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['Galaxy with f_esc, HI=1, HeI=1'] = fesc1

    filename_fesc = cwd+'data/fesc/fe_lyccontrot0.001hi1.000hei.txt'
    tab = ascii.read(filename_fesc)
    sp = S.ArraySpectrum(wave=tab['lam'], flux=tab['lh1=17.5'], waveunits='Angstrom', fluxunits='flam')
    fesc1 = sp.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['Galaxy with f_esc, HI=0.001, HeI=1'] = fesc1

    flatsp = S.FlatSpectrum(21, fluxunits='flam')
    flat = flatsp.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    spec_dict['Flat in F_lambda'] = flat  


    return spec_dict

