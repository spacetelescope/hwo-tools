;  stellar flux function
;    lam - wavelength (um)
;   Teff - effective temperature (K)
;     Rs - stellar radius (solar radii)
;      d - distance to star (pc)
;     AU - flag that indicates d is in AU
;  Fstar - stellar flux (W/m**2/um)
FUNCTION Fstar, lam, Teff, Rs, d, AU=au
  Rsun  = 6.958d8       ; solar radius (m)
  ds    = 3.08567d16    ; parsec (m)
  IF KEYWORD_SET( AU ) THEN ds = 1.495979d11 ; AU (m)
  lambda= 1.d-6*lam     ; wavelength (m)
  c1    = 3.7417715d-16 ; 2*pi*h*c*c (kg m**4 / s**3)   
  c2    = 1.4387769d-2  ; h*c/k (m K)
  pow   = c2/lambda/Teff
  Fs    = c1/( (lambda^5.)*(exp(pow)-1.d) )*1.d-6
  RETURN, Fs*(Rs*Rsun/d/ds)^2.
END

;  planetary flux function
;      A - planetary geometric albedo
;    Phi - planetary phase function
;  Fstar - stellar flux (W/m**2/um)
;     Rp - planetary radius (Earth radii)
;      d - distance (pc)
;     au - flag that indicates d is in AU
;  Fplan - planetary flux (W/**2/um)
FUNCTION Fplan, A, Phi, Fstar, Rp, d, AU=au
  Re    = 6.371e6       ; radius of Earth (m)
  ds    = 3.08567d16    ; parsec (m)
  IF KEYWORD_SET( AU ) THEN ds = 1.495979d11 ; AU (m)
  RETURN, A*Phi*Fstar*(Rp*Re/d/ds)^2.
END

;  planet-star flux ratio
;      A - planetary geometric albedo
;    Phi - planetary phase function
;     Rp - planetary radius (Earth radii)
;      r - orbital distance (au)
;   FpFs - planet-star flux ratio
FUNCTION FpFs, A, Phi, Rp, r
  Re = 6.371e6 ; radius of Earth (m)
  ds = 1.495979d11 ; AU (m)
  RETURN, A*Phi*(Rp*Re/r/ds)^2.
end

;  planet photon count rate
;      q - quantum efficiency
;    fpa - fraction of planet flux in Airy pattern
;      T - system throughput
;    lam - wavelength (um)
;   dlam - spectral element width (um)
;  Fplan - planetary flux (W/m**2/um)
;      D - telescope diameter (m)
;  cplan - planet photon count rate (s**-1)
FUNCTION cplan, q, fpa, T, lam, dlam, Fplan, D
  hc  = 1.986446d-25 ; h*c (kg*m**3/s**2)
  RETURN, !DPI*q*fpa*T*(lam*1.d-6/hc)*dlam*Fplan*(D/2)^2.
END

;  zodiacal light count rate
;      q - quantum efficiency
;      X - diameter or length of photometric aperture (lambda/D)
;      T - system throughput
;    lam - wavelength (um)
;   dlam - spectral element width (um)
;      D - telescope diameter (m)
;    MzV - zodiacal light surface brightness (mag/arcsec**2)
;   /SUN - set to use WMO solar spectrum
;  /CIRC - keyword to use a circular aperture
;  czodi - zodiacal light photon count rate (s**-1)
FUNCTION czodi, q, X, T, lam, dlam, D, MzV, SUN=sun, CIRC=circ
  hc    = 1.986446d-25 ; h*c (kg*m**3/s**2)
  F0V   = 3.6d-8       ; zero-mag V-band flux (W/m**2/um)
  FsolV = 1.86d+3      ; Solar V-band flux at 1 AU
  IF KEYWORD_SET(sun) THEN BEGIN
    fn    = './input_data/wmo_solar_spectrum.dat'
    READCOL, fn, lamsol, Fsol, SKIPLINE=32, /SILENT
    Fsol  = DEGRADE_SPEC(Fsol,lamsol,lam,DLAM=dlam) ; degrade solar spectrum (W/m**2/um)
  ENDIF ELSE BEGIN
   Teffs  = 5778.d ;Sun effective temperature
      Rs  = 1.d    ;Sun radius (in solar radii)
    Fsol  = Fstar(lam, Teffs, Rs, 1., /AU)
  ENDELSE
  rat   = DBLARR(N_ELEMENTS(lam))
  rat[*]= Fsol[*]/FsolV ; ratio of solar flux to V-band solar flux
  IF KEYWORD_SET(circ)  THEN Omega = !DPI*(X/2.*lam*1.d-6/D*180.*3600./!DPI)^2. ; circular aperture size (arcsec**2)
  IF ~KEYWORD_SET(circ) THEN Omega = (X*lam*1.d-6/D*180.*3600./!DPI)^2. ; square aperture size (arcsec**2)
  RETURN, !DPI*q*T*Omega*dlam*(lam*1.d-6/hc)*(D/2)^2.*rat*F0V*10^(-Mzv/2.5)
END

;  exozodiacal light count rate
;      q - quantum efficiency
;      X - diameter or length of photometric aperture (lambda/D)
;      T - system throughput
;    lam - wavelength (um)
;   dlam - spectral element width (um)
;      D - telescope diameter (m)
;      r - orbital distance (au)
;  Fstar - host star spectrum *at 1 au* (W/m**2/um)
;    Nez - number of exozodis
;   MezV - exozodiacal light surface brightness (mag/arcsec**2)
;   /SUN - set to use WMO solar spectrum
;  /CIRC - keyword to use a circular aperture
; cezodi - exozodiacal light photon count rate (s**-1)
FUNCTION cezodi, q, X, T, lam, dlam, D, r, Fstar, Nez, MezV, SUN=sun, CIRC=circ
  hc    = 1.986446d-25 ; h*c (kg*m**3/s**2)
  F0V   = 3.6d-8       ; zero-mag V-band flux (W/m**2/um)
  FsolV = 1.86d+3      ; Solar V-band flux at 1 AU
  IF KEYWORD_SET(sun) THEN BEGIN
    fn    = './input_data/wmo_solar_spectrum.dat'
    READCOL, fn, lamsol, Fsol, SKIPLINE=32, /SILENT
    Fsol  = DEGRADE_SPEC(Fsol,lamsol,lam,DLAM=dlam) ; degrade solar spectrum (W/m**2/um)
  ENDIF ELSE BEGIN
   Teffs  = 5778.d ;Sun effective temperature
      Rs  = 1.d    ;Sun radius (in solar radii)
    Fsol  = Fstar(lam, Teffs, Rs, 1., /AU) ; Sun as blackbody (W/m**2/um)
  ENDELSE
  rat   = DBLARR(N_ELEMENTS(lam))
  rat[*]= Fstar[*]/FsolV ; ratio of solar flux to V-band solar flux
  IF KEYWORD_SET(circ)  THEN Omega = !DPI*(X/2.*lam*1.d-6/D*180.*3600./!DPI)^2. ; circular aperture size (arcsec**2)
  IF ~KEYWORD_SET(circ) THEN Omega = (X*lam*1.d-6/D*180.*3600./!DPI)^2. ; square aperture size (arcsec**2)
  RETURN, !DPI*q*T*Omega*dlam*(lam*1.d-6/hc)*(D/2)^2.*(1./r)^2.*rat*Nez*F0V*10^(-Mezv/2.5)
END

;  speckle count rate
;      q - quantum efficiency
;      T - system throughput
;      C - design contrast
;    lam - wavelength (um)
;   dlam - spectral element width (um)
;      D - telescope diameter (m)
;  Fstar - host star spectrum at distance to system (W/m**2/um)
; cspeck - speckle photon count rate (s**-1)
FUNCTION cspeck, q, T, C, lam, dlam, Fstar, D
  hc    = 1.986446d-25 ; h*c (kg*m**3/s**2)
  RETURN, !DPI*q*T*C*dlam*Fstar*(lam*1.d-6/hc)*(D/2)^2.
END

;  dark count rate
;     De - dark count rate (s**-1)
;      X - diameter or length of photometric aperture (lambda/D)
;    lam - wavelength (um)
;      D - telescope diameter (m)
;  theta - angular diameter of lenslet or pixel (arcsec)
; DNhpix - number of pixels spectrum spread over in horizontal, for IFS
; /IMAGE - keyword set to indicate imaging mode (not IFS)
;  /CIRC - keyword to use a circular aperture
;  cdark - dark count rate (s**-1)
FUNCTION cdark, De, X, lam, D, theta, DNhpix, IMAGE=image, CIRC=circ
  IF KEYWORD_SET(circ)  THEN Omega = !DPI*(X*lam*1.d-6/D*180.*3600./!DPI)^2. ; circular aperture diameter (arcsec**2)
  IF ~KEYWORD_SET(circ) THEN Omega = 4.d0*(X*lam*1.d-6/D*180.*3600./!DPI)^2. ; square aperture diameter (arcsec**2)
  Npix  = Omega/!DPI/theta^2.
  IF ~KEYWORD_SET( IMAGE ) THEN Npix = 2*DNhpix*Npix
  RETURN, De*Npix
END

;  read noise count rate
;     Re - read noise counts per pixel
;      X - diameter or length of photometric aperture (lambda/D)
;    lam - wavelength (um)
;      D - telescope diameter (m)
;  theta - angular diameter of lenslet or pixel (arcsec)
;  Dtmax - maximum exposure time (hr)
; /IMAGE - keyword set to indicate imaging mode (not IFS)
;  /CIRC - keyword to use a circular aperture
;  cread - read count rate (s**-1)
FUNCTION cread, Re, X, lam, D, theta, DNhpix, Dtmax, IMAGE=image, CIRC=circ
  IF KEYWORD_SET(circ)  THEN Omega = !DPI*(X*lam*1.d-6/D*180.*3600./!DPI)^2. ; circular aperture diameter (arcsec**2)
  IF ~KEYWORD_SET(circ) THEN Omega = 4.d0*(X*lam*1.d-6/D*180.*3600./!DPI)^2. ; square aperture diameter (arcsec**2)
  Npix  = Omega/!DPI/theta^2.
  IF ~KEYWORD_SET( IMAGE ) THEN Npix = 2*DNhpix*Npix
  RETURN, Npix/(Dtmax*3600.)*Re
END

;  clock induced charge count rate
;     Rc - clock induced charge counts per pixel per read
;      X - diameter or length of photometric aperture (lambda/D)
;    lam - wavelength (um)
;      D - telescope diameter (m)
;  theta - angular diameter of lenslet or pixel (arcsec)
;  Dtmax - maximum exposure time (hr)
; /IMAGE - keyword set to indicate imaging mode (not IFS)
;  /CIRC - keyword to use a circular aperture
;  cread - read count rate (s**-1)
FUNCTION ccic, Rc, X, lam, D, theta, DNhpix, Dtmax, IMAGE=image, CIRC=circ
  IF KEYWORD_SET(circ)  THEN Omega = !DPI*(X*lam*1.d-6/D*180.*3600./!DPI)^2. ; circular aperture diameter (arcsec**2)
  IF ~KEYWORD_SET(circ) THEN Omega = 4.d0*(X*lam*1.d-6/D*180.*3600./!DPI)^2. ; square aperture diameter (arcsec**2)
  Npix  = Omega/!DPI/theta^2.
  IF ~KEYWORD_SET( IMAGE ) THEN Npix = 2*DNhpix*Npix
  RETURN, Npix/(Dtmax*3600.)*Rc
END

;  fraction of Airy power contained in square or circular aperture
;      X - size of photometric aperture (lambda/D)
;  /CIRC - keyword to use a circular aperture
; f_airy - fraction of power in Airy pattern of size X*lambda/D
FUNCTION f_airy, X, CIRC=circ
  IF KEYWORD_SET(circ) THEN BEGIN
    fpa = 1. - BESELJ(!DPI*X,0)^2. - BESELJ(!DPI*X,1)^2. ;fraction of power in circular Airy pattern
  ENDIF ELSE BEGIN
    X_grid   = DINDGEN(100)/10.d
    fpa_grid = $                 ;fpa on grid, computed using code below
    [ 0.00000,    0.0309038,     0.117701,     0.244285,     0.388709,     0.528808, $
     0.647304,     0.734911,     0.790813,     0.820765,     0.833866,     0.839249, $
     0.843710,     0.850753,     0.860930,     0.872944,     0.884901,     0.895255, $
     0.903243,     0.908880,     0.912712,     0.915504,     0.917992,     0.920720, $
     0.923959,     0.927698,     0.931686,     0.935540,     0.938892,     0.941519, $
     0.943427,     0.944836,     0.946081,     0.947471,     0.949178,     0.951186, $
     0.953330,     0.955383,     0.957163,     0.958591,     0.959709,     0.960639, $
     0.961530,     0.962501,     0.963607,     0.964837,     0.966128,     0.967391, $
     0.968543,     0.969532,     0.970352,     0.971048,     0.971696,     0.972377, $
     0.973147,     0.974014,     0.974941,     0.975857,     0.976696,     0.977418, $
     0.978026,     0.978559,     0.979075,     0.979627,     0.980241,     0.980912, $
     0.981611,     0.982297,     0.982935,     0.983505,     0.984008,     0.984464, $
     0.984904,     0.985360,     0.985853,     0.986387,     0.986949,     0.987510, $
     0.988040,     0.988519,     0.988944,     0.989331,     0.989707,     0.990099, $
     0.990524,     0.990983,     0.991459,     0.991931,     0.992376,     0.992783, $
     0.993154,     0.993502,     0.993844,     0.994197,     0.994571,     0.994966, $
     0.995376,     0.995784,     0.996176,     0.996542]
    ;
    fpa = INTERPOL(fpa_grid,X_grid,X)
  ENDELSE
  RETURN, fpa
END

; numerical integration to determine fraction of Airy power contained square aperture (SLOW!)
;      X - size of photometric aperture (lambda/D)
; f_airy - fraction of Airy power in aperture of size X*lambda/D
FUNCTION f_airy_int, X
  N     = 1000 ; sets horizontal and vertical grid resolution
  E0    = 1.27324 ; total power contained in Airy pattern
  xg = X*FINDGEN(N)/(N-1) ; horizontal grid
  yg = X*FINDGEN(N)/(N-1) ; vertical grid
  Iairy = DBLARR(N,N)
  FOR i=0,N-1 DO FOR j=0,N-1 DO $
  Iairy[i,j] = 4*BESELJ(!DPI*SQRT(xg[i]^2. + yg[j]^2.),1)^2./ $
               (!DPI*SQRT(xg[i]^2. + yg[j]^2.))^2. ; compute Airy intensity
  Iairy[0,0] = 1.d0
  E     = 0.d0
  dA    = (X/N)^2. ; area of each dxdy element
  FOR i=0,N-1 DO FOR j=0,N-1 DO $
    E = E + Iairy[i,j]*dA ; sum intensities
  E     = 4.d*E ; factor of 4 as integral only over one quadrant
  fpa   = E/E0  
END