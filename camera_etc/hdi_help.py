def help(): 
   return """
      <div class="container"> 
        <div class="col-lg-6">
         <p align="justify">This is the basic ETC for photometry in multiband images. 
            Choose your telescope aperture, exposure time, and magnitude normalization. 
            The normalization is done in the V band (550 nm). 
        <p> 
        <p> 
        <p align="justify"> Given an aperture and magnitude, choose the exposure time that reaches your desired S/N.  
        <p>  
        <p>  
        <p align="justify">To obtain limiting magnitudes given exposure time, set that time and then tune the magnitude to reach your desired limiting S/N.</p>
        </div>
        <div class="col-lg-6">
        <p>  
        <p>  
          <p align="justify">The details of these calculations are <a href="http://jt-astro.science/luvoir_simtools/hdi_etc/SNR_equation.pdf" target="_blank"> here</a>. 
             We assume that the pixel size in each band critically samples the telescope's diffraction limited PSF at the shortest wavelength in that channel. 
             Thermal backgrounds are computed for T = 280 K, which substantially affects the K band.</p>
        <p>
          <p align="justify">Uploaded spectra can be in either fixed width two-column ASCII (wave flux) or FITS format, where the spectrum is in the second HDU (HDU1, BINTABLE) with column 1 = "WAVELENGTH" and column 2 = "FLUX". Wavelength is assumed to be in Angstroms, Flux in FLAM (erg s−1 cm−2˚𝐴−1)</p>
        </div>
      </div>
          """

