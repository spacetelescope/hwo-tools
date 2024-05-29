

def uvspec_snr(telescope, mode, template, fuvmag, exptime): 
    ''' Run a basic SNR calculation that takes in a telescope, 
        spectral template, normalization magnitude, and exposure 
        time to compute SNR
        This function is supplied as a user-friendly wrapper to 
        the SYOTools API '''  
    
    from syotools.models import Telescope, Spectrograph
    from syotools.utils.jsonunit import str_jsunit
    from syotools.spectra import SpectralLibrary
    import astropy.units as u

    # create the basic objects 
    uvi, tel = Spectrograph(), Telescope() 
    tel.set_from_json(telescope)
    tel.add_spectrograph(uvi)
    uvi_exp = uvi.create_exposure()
    uvi.mode = mode

    uvi_exp.sed_id = template
    uvi_exp.renorm_sed(fuvmag * u.ABmag, bandpass='galex,fuv')

    #also set the exposure time to 5 hours for all modes 
    # JT ------- IT IS NOT CLEAR WHAT THESE ELEMENTS CORRESPOND TO 
    uvi_exp.exptime[1]['value'] = [exptime, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]

    #Print the current template & mode
    print("Current SED template: {}".format(SpectralLibrary[uvi_exp.sed_id])) 
    print("Current grating mode: {}".format(uvi.descriptions[uvi.mode])) 
    print("Current exposure time: {} hours\n".format(uvi_exp.exptime[1]['value'][0])) 
    
    lumos_template_codes = ['flam', 'qso', 's99', 'o5v', 'g2v', 'g191b2b', 'gd71', 'gd153', 'ctts', 
                        'mdwarf', 'orion', 'nodust', 'ebv6', 'hi1hei1', 'hi0hei1']
    available_lumos_templates = [SpectralLibrary[tc] for tc in lumos_template_codes]

    snr = uvi_exp.recover('snr')
    uvi_exp.enable()
    uvi_sed, uvi_snr = uvi_exp.recover('sed', 'snr')

    wave, snr =  uvi.wave[1]['value'], uvi_exp.snr[1]['value']

    return wave, snr, uvi 