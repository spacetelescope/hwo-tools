
def uvspec_snr(telescope, mode, template, fuvmag, exptime, silent=True): 
    ''' Run a basic SNR calculation that takes in a telescope, 
        spectral template, normalization magnitude, and exposure 
        time to compute SNR. For converting magnitude, template, 
	      and SNR to a desired exposure time, use uvspec_exptime.py 
      
        usage: 
	      wave, snr, uvi = uvspec_snr(telescope, mode, template, uvmag, exptime) 

          positional arguments:

           1-telescope = 'EAC1', 'EAC2', or 'EAC3'. This argument is a string. 
             EAC1 = 6 m inner diameter, 7.2 outer diameter hex pattern, off-axis 
             EAC2 = 6 m diameter off-axis 
             EAC3 = 8 m diameter on-axis 

           2-mode = your choice of UVI grating, a string: 
		        ['G120M', 'G150M', 'G180M', 'G155L', 'G145LL', 'G300M']

           3-template = your choice of spectral template: 
		          ['flam', 'qso', 's99', 'o5v', 'g2v', 'g191b2b', 'gd71', 'gd153', 'ctts', 
                        'mdwarf', 'orion', 'nodust', 'ebv6', 'hi1hei1', 'hi0hei1']

           4-fuvmag = FUV magnitude to normalize the template spectrum, a float. 

	   5-exptime = desired exposure time in hours, a float 

        outputs are two arrays of floats for wavelength and snr and the Spectrograph 
		    object in case it is needed by other code. 
       '''  
    
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

    uvi_exp.exptime[1]['value'] = [exptime, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]

    uvi_exp.verbose = True 
    tel.verbose = True 
    if (silent):  
       uvi_exp.verbose = False 
       tel.verbose = False 
       uvi.verbose = False 
       print("We have set verbose = False") 

    if not silent: 
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

def uvspec_exptime(telescope, mode, template, fuvmag, snr_goal, silent=True): 
    ''' Run a basic SNR calculation that takes in a telescope, 
      spectral template, normalization magnitude, and SNR goal  
      to compute exposure time. For converting magnitude, template, 
      and exptime to SNR, use uvspec_snr.py 
      
        usage: 
	       wave, exptime, uvi = uvspec_exptime(telescope, mode, template, uvmag, snr_goal) 

          positional arguments:

           1-telescope = 'EAC1', 'EAC2', or 'EAC3'. This argument is a string. 
             EAC1 = 6 m inner diameter, 7.2 outer diameter hex pattern, off-axis 
             EAC2 = 6 m diameter off-axis 
             EAC3 = 8 m diameter on-axis 

           2-mode = your choice of UVI grating, a string: 
		        ['G120M', 'G150M', 'G180M', 'G155L', 'G145LL', 'G300M']

           3-template = your choice of spectral template: 
	          	['flam', 'qso', 's99', 'o5v', 'g2v', 'g191b2b', 'gd71', 'gd153', 'ctts', 
                        'mdwarf', 'orion', 'nodust', 'ebv6', 'hi1hei1', 'hi0hei1']

           4-fuvmag = FUV magnitude to normalize the template spectrum, a float. 

	         5-snr_goal = desired SNR, per pixel 

         outputs are two arrays of floats for wavelength and exptime and the Spectrograph 
		     object in case it is needed by other code. 
       '''  
    
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

    uvi_exp.verbose = True 
    if (silent):  uvi_exp.verbose = False 

    if not silent: 
        print("Current SED template: {}".format(SpectralLibrary[uvi_exp.sed_id])) 
        print("Current grating mode: {}".format(uvi.descriptions[uvi.mode])) 
        print("Current exposure time: {} hours\n".format(uvi_exp.exptime[1]['value'][0])) 
    
    lumos_template_codes = ['flam', 'qso', 's99', 'o5v', 'g2v', 'g191b2b', 'gd71', 'gd153', 'ctts', 
                        'mdwarf', 'orion', 'nodust', 'ebv6', 'hi1hei1', 'hi0hei1']
    available_lumos_templates = [SpectralLibrary[tc] for tc in lumos_template_codes]

    uvi_exp._snr_goal= snr_goal * (u.ct)**0.5 / (u.pix)**0.5

    snr = uvi_exp.recover('exptime')
    uvi_exp.unknown = 'exptime' #< --- this triggers the _update_exptime function in the SpectrographicExposure exposure object 

    #uvi_exp.enable()
    uvi_sed, uvi_exptime = uvi_exp.recover('sed', 'exptime')

    wave, exptime =  uvi.wave[1]['value'], uvi_exp.exptime[1]['value']

    return wave, exptime, uvi 
