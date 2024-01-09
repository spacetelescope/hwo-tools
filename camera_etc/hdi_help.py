

def help(): 
   return """

      <div class="container"> 
        <div class="col-lg-6">
         <p>This is the basic ETC for photometry in multiband images. Choose your telescope aperture, exposure time, and magnitude normalization. The normalization is done in the V band (550 nm). 
        <p> 
        <p> Given an aperture and magnitude, choose the exposure time that reaches your desired S/N.  
        <p> 
        <p>To obtain limiting magnitudes given exposure time, set that time and then tune the magnitude to reach your desired limiting S/N.</p>
        </div>
        <div class="col-lg-6">
        <p> 
          <p>The details of these calculations are <a href="http://jt-astro.science/luvoir_simtools/hdi_etc/SNR_equation.pdf" target="_blank"> here</a>. We assume that the pixel size in each band critically samples the telescope's diffraction limited PSF at the shortest wavelength in that channel. Thermal backgrounds are computed for T = 280 K, which substantially affects the K band.</p>
        </div>
      </div>
          """

