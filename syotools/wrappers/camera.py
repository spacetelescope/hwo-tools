def camera_snr(telescope, template, mag, exptime): 
    
    from syotools.models import Telescope, Camera
    from syotools.utils.jsonunit import str_jsunit
    from syotools.spectra import SpectralLibrary
    import astropy.units as u
    import numpy as np 

    # create the basic objects 
    hri, tel = Camera(), Telescope() 
    tel.set_from_json(telescope)
    tel.add_camera(hri)
    hri_exp = hri.create_exposure()

    hri_exp.sed_id = template
    hri_exp.renorm_sed(mag * u.ABmag, bandpass='v')
    exptime_list = []
    for i in hri.bandnames: exptime_list.append(exptime) 
    hri_exp.exptime[1]['value'] = exptime_list 

    #Print the current template & mode
    print("Current SED template: {}".format(hri_exp.sed_id)) 
    print("Current exposure time: {} hours\n".format(hri_exp.exptime[1]['value'][0])) 
    
    snr = hri_exp.recover('snr')
    hri_exp.enable()
    hri_sed, hri_snr = hri_exp.recover('sed', 'snr')  
    snr = hri_exp.snr[1]['value']

    for bb, ss in zip(hri.bandnames, snr): print("{}, SNR = {}".format(bb, ss)) 

    return snr, hri
