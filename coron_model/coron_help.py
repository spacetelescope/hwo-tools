

def help(): 
   return """

      <div class="container">
            <p>The "Planet" tab includes options to simulate several types of planetary spectra that can be selected from the "Planet Spectrum" dropdown menu. The telescope-planetary system separation distance can be set using the "Distance" slider. When a target is selected, the "Planet Radius" and "Semi-major axis of orbit" sliders will default to the correct positions for the selected planet. Note that while it is possible to adjust these parameters for each target, changing them can result in spectra representing non-physical targets. Also included under the "Planet" tab is a slider for scaling exozodiacal dust. 
              
              <p>The "Observation" tab controls telescope integration time per coronagraphic bandpass, mirror diameter, spectrograph resolution for UV-VIS-NIR channels, telescope temperature, and the option to turn on a ground-based simulator.
              
              <p>The "Instrumentation" tab controls the instrument inner working angle (IWA), outer working angle (OWA), both in terms of lambda/D, and the maximum time for a single exposure.
              
              <p>The "Exposure Time Calculator" tab contains a slider to set a desired signal-to-noise ratio. In the "Exposure Time" plot tab, the simulator will display the integration time required to obtain this signal-to-noise ratio for the current telescope and instrumentation setup. Note that this tab applies only to the Exposure Time plot, not to the Spectrum plot. 

              <p>In the "Download" tab, spectral data can be downloaded in either .txt or .fits format.

              <p>The underlying model is derived from the python-based version of Tyler Robinson's coronagraphic spectrum and noise model (Robinson et al. 2016). Python by Jacob Lustig-Yaeger. Bokeh rendering by Jason Tumlinson and Giada Arney. </p>
              For full details, please see the readme file <a href="http://jt-astro.science/luvoir_simtools/coron_model/coron_readme.txt" target="_blank">here</a>.
      </div>


          """

