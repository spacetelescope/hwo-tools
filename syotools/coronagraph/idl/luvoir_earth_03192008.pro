PRO LUVOIR_EARTH_03192008

  ;;; set integration time
  ;
  Dthr     = 200.    ;integration time (hr)
  Dt     = Dthr*3600 ;convert to (s)

  ;;; set planet parameters
  ;
  ;   spectrum
  tag = 'earth_03192008'
  savdir = '/astro/users/jlustigy/Dropbox/Astronomy/Plots/coronagraph/data_files/'
  fn    = 'spectra/earth_03192008_px34_hitran2012_50_100000cm_toa.rad'
  READCOL, fn, lamhr, wnhr, solhr, flxhr, r1, r2, r3, r4   ;<- read SMART output
  refl  = flxhr/solhr ;reflectance
  Ahr   = GAUSSFOLD(ALOG(lamhr),refl,1./70.) ;Gaussian convolution to R=70

  ;   other planet params
  alpha = 90.     ;phase angle at quadrature
  Phi   = 1.      ;phase function at quadrature (already included in SMART run)
  Rp    = 1.0     ;Earth radii
  r     = 1.0     ;semi-major axis (AU)
  
  ;;; set star parameters
  ;
  Teff  = 5780.   ;Sun-like Teff (K)
  Rs    = 1.      ;star radius in solar radii
  
  ;;; set planetary system parameters
  ;
  d    = 10.     ;distance to system (pc)
  Nez  = 1.      ;number of exo-zodis

  ;;; call noise model
  ;
  CORONAGRAPH_LUVOIR, Ahr, lamhr, solhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez, $
                 lammin, lammax, Res, X, fpa, diam, Tput, C, IWA, OWA, $
                 lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, $
                 ctot, DtSNR, /COMPUTE_LAM, /COMPUTE_DLAM, /NIR

  ;;; compute SNR and 1-sigma
  ;
  SNR  = cp*Dt/SQRT((cp + 2*(csp+cz+cez+cD+cR))*Dt)
  sigma= Cratio/SNR

  ;;; make noised-up spectrum
  ;
  cont = Cratio + RANDOMN(seed,N_ELEMENTS(lam))*sigma ; noisy planet-star contrast
  i    = WHERE(cont LT 0)
  j    = WHERE(sigma EQ !VALUES.D_INFINITY)
  i2   = WHERE(cont EQ !VALUES.D_INFINITY)
  IF(i[0] NE -1) THEN cont[i] = 0.0 ;<- spectrum cannot fall below zero
  IF(j[0] NE -1) THEN BEGIN 
    sigma[j] = !VALUES.F_NAN ;<- Infinite Errors!
    cont[j] = !VALUES.F_NAN
  ENDIF
  IF(i2[0] NE -1) THEN cont[i2] = !VALUES.F_NAN ;<- No planet counts 
  IF(i2[0] NE -1) THEN sigma[i2] = !VALUES.F_NAN ;<- Infinite Errors! 

  ;;; make a plot
  ;
  ;   some sizing parameters
  cth = 5.0   ;character thickness
  csi = 0.9  ;character size
  th  = 6.    ;line thickness
  colorFile = 'fsc_brewer.tbl'
  v = FINDGEN(17) * (!PI*2/16.)
  USERSYM, COS(v), SIN(v), /FILL

  ;
  SET_PLOT, 'ps'
  DEVICE, FILENAME='plots/luvoir_'+tag+'.pdf'
  DEVICE, /COLOR
  LOADCT, 0
  PLOT, lam, Cratio*1.e9, XRANGE=[lammin,2.0], XSTYLE=1, $
        YRANGE=[0.,1.05*MAX(Cratio)*1.e9], $
        YSTYLE=1, YTHICK=th, CHARSIZE=csi, XTITLE=TeXtoIDL('Wavelength (\mum)'), $
        YTITLE=TeXtoIDL('Planet-Star Flux Ratio \times 10^{9}'), XTHICK=th, $
        POSITION=[0.155,0.195,0.95,0.895], CHARTHICK=cth, /NODATA;, /XLOG
  LOADCT, 15
  OPLOT, lam, Cratio*1.e9, THICK=th, COLOR=40, PSYM=10
  LOADCT, 0
  OPLOTERROR, lam, cont*1.e9, sigma*1.e9, THICK=0.75*th, ERRCOLOR=125, PSYM=3
  LOADCT, 15
  OPLOT, lam, cont*1.e9, PSYM=8, SYMSIZE=0.1*th, COLOR=175
  DEVICE, /CLOSE
  SET_PLOT, 'X'
  LOADCT, 0

  ;;; save data file: wavelength, true planet/star photon counts, noised planet/star photon counts, 1-sigma errors
  WRITECOL, savdir+'luvoir_'+tag+'.txt', lam, Cratio, cont, sigma

END
