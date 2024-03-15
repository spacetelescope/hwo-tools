;+
; NAME:     
;           gaussfold
;
;
; PURPOSE:
;           Smoothes a plot by convolving with a Gaussian profile.
;           Main purpose is to convolve a spectrum (flux against
;           wavelength) with a given instrument resolution.
;           Also applicable e.g. to smooth ligthcurves or in
;           time-series analysis. 
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;           smoothedFlux = gaussfold(lambda, flux, sigma,
;                                    LAMMIN=lammin, LAMMAX=lammax) 
;
;
; INPUTS:
;           lambda = In ascending order sorted array (double or float)
;                    containing the values of the x-axis (e.g. the
;                    wavelength- or frequencygrid). Irregularly spaced
;                    grids are acepted. 
;           flux   = Array (double or float) same size as lambda
;                    containing the values of the y-axis (e.g. the
;                    flux).  
;           fwhm   = FWHM in units of the x-axis (e.g. Angstroem) of
;                    the Gaussian profile.
;
;
; OPTIONAL INPUTS:
;           NONE
;
;
; KEYWORD PARAMETERS:
;           lammin &
;           lammax   = Defines the x-axis range. All y-data within
;                      this range will be smoothed.  
;                      CAUTION: improtant in case of large arrays (memory!)
;                      DEFAULTS: lammin = MIN(lambda) i.e. lambda[0]
;                                lammax = MAX(lambda) i.e. lambda[N_ELEMENTS(lambda)-1]
;
;
; OUTPUTS:
;           smoothedFlux = Array (double of float) same size as lambda
;                          containing the smoothed flux.  
;
;
; OPTIONAL OUTPUTS:
;           NONE
;
;
; COMMON BLOCKS:
;           NONE
;
;
; SIDE EFFECTS:
;           NONE
;
;
; RESTRICTIONS:
;           NONE
;
;
; PROCEDURE:
;           1. Interpolation of the flux on a fine spaced grid.
;              Oversamplingfaktor: 17
;           2. Convolution of the flux with a gaussian profile 
;           3. Interpolation of the smoothed flux on the original grid.
;
;
; EXAMPLE:
;           fluxS = GAUSSFOLD( lambda, flux, 0.5)
;
;
; MODIFICATION HISTORY:
;           Version 1.0  by Katja Pottschmidt  1999
;           Version 2.0  by Jochen Deetjen     1999
;                        - Documentation added
;                        - Keyword defaults defined
;                        - convolution with the gaussian kernel is
;                          restricted within a window 
;                          -> by far less time consuming 
;
;-


FUNCTION gaussfold, lam, flux, fwhm, lammin=lammin, lammax=lammax
   
   IF (NOT KEYWORD_SET(lammin)) THEN $
     lammin = MIN(lam)
   
   IF (NOT KEYWORD_SET(lammax)) THEN $
     lammax = MAX(lam)
   
   dlambda   = fwhm / 17D0
   
   interlam  = lammin + dlambda * DINDGEN( LONG( (lammax-lammin)/dlambda+1 ))
   interflux = INTERPOL( flux, lam, interlam )   
   
   fwhm_pix = fwhm / dlambda
   window   = FIX( 17 * fwhm_pix )
   
   ;; get a 1D gaussian profile
   gauss = PSF_GAUSSIAN( NP=window, FWHM=fwhm_pix, /NORMALIZE, NDIMEN=1 )
   
   ;; convolve input spectrum with the gauss profile
   fold  = CONVOL( interflux, gauss, /CENTER, /EDGE_TRUNCATE, /NAN)

   fluxfold  = INTERPOL( fold, interlam, lam )
   
   RETURN, fluxfold
   
END










