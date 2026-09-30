import syotools

def help(): 
   return f"""

      <div class="container"> 
          <p align="justify">This is the basic ETC for UV spectroscopy with HWO's UV spectrograph, powered by SYOTools v{syotools.__version__}.<br> 
          <p align="justify">The top controls model the source spectrum with various templates,  
                an optional redshift, and magnitude normalization. All template spectra are normalized to the given magnitude in the GALEX FUV band. 
          <p align="justify">The top plot shows the input source spectrum in dark red  
                and the "Background Equivalent Flux" for that mode. 
          <p align="justify">The lower panel shows the signal-to-noise ratio for the selected grating, aperture, and exposure time. 
        <p>
          <p align="justify">Uploaded spectra can be in either fixed width two-column ASCII (wave flux) or FITS format, where the spectrum is in the second HDU (HDU1, BINTABLE) with column 1 = "WAVELENGTH" and column 2 = "FLUX". Wavelength is assumed to be in Angstroms, Flux in FLAM (erg s−1 cm−2˚𝐴−1)</p>

      </div>
          """

