;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Quick LUVOIR noise model - Tyler D. Robinson
;  inputs:
;    Ahr  - hi-res planetary albedo spectrum
;  lamhr  - wavelength grid for Ahr (um)
;  solhr  - hi-res TOA solar spectrum (W/m**2/um)
;  alpha  - phase angle (deg)
;    Phi  - phase function evaluated at alpha
;     Rp  - planetary radius (Rearth)
;   Teff  - stellar effective temperature (K)
;     Rs  - stellar radius (Rsun)
;      r  - orbital separation (au)
;      d  - distance to system (pc)
;    Nez  - number of exozodis
;  lammin - minimum wavelength (um)
;  lammax - maximum wavelength (um)
;     Res - spectral resolution (lambda/Dlambda)
;    diam - telescope diameter (m)
;    Tput - system throughput
;      C  - raw contrast
;     IWA - inner working angle (lambda/D)
;     OWA - outer working angle (lambda/D; unless /FIXOWA)
;
;  outputs:
;    lam  - low-res wavelength grid (um)
;   dlam  - spectral bin widths (um)
;      A  - planetary albedo spectrum at low-res
;      q  - quantum efficiency
; Cratio  - planet-star contrast (flux) ratio
;     cp  - planet count rate (s**-1)
;    csp  - speckle count rate (s**-1)
;     cz  - zodiacal light count rate (s**-1)
;    cez  - exozodiacal light count rate (s**-1)
;     cD  - dark current count rate (s**-1)
;     cR  - read noise count rate (s**-1)
;   ctot  - total count rate (s**-1)
;  DtSNR  - integration time to SNR=1 (hr)
;
;  options:
;        FIX_OWA - set to fix OWA at OWA*lammin/D, as would occur
;                  if lenslet array is limiting the OWA
;    COMPUTE_LAM - set to compute lo-res wavelength grid, otherwise
;                  the grid input as variable 'lam' is used
;   COMPUTE_DLAM - set to compute lo-res wavelength grid widths,
;                  otherwise the grid input as variable 'dlam' is used
;            NIR - re-adjusts pixel size in NIR, as would occur if a
;                  second instrument was designed to handle the NIR
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO CORONAGRAPH_LUVOIR, Ahr, lamhr, solhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez, $
                 lammin, lammax, Res, X, fpa, diam, Tput, C, IWA, OWA, $
                 lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, $
                 ctot, DtSNR, FIX_OWA = fix_owa, COMPUTE_LAM = compute_lam, $
                 COMPUTE_DLAM = compute_dlam, NIR = nir

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; set key LUVOIR parameters ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  diam   = 10.     ;diameter (m)
  C      = 1.d-10  ;raw contrast
  lammin = 0.4     ;minimum wavelength (um)
  lammax = 2.5     ;maximum wavelength (um)
  Res    = 70.     ;spectral resolution (lambda/Dlambda)
  Tput   = 0.05    ;system throughput
  IWA    = 3.      ;inner working angle (lambda/D)
  OWA    = 20.     ;outer working angle (lambda/D)
  Tsys   = 150.    ;telescope/system temperature (K)
  Tdet   = 50.     ;detector temperature (K)
  emis   = 0.9     ;telescope/system emissivity
  X      = 1.5     ;aperture size (lambda/D)
  fpa    = f_airy(X) ; fraction of planetary signal in Airy pattern
  DNhpix = 3      ; horizontal pixel spread of IFS spectrum
  Dtmax  = 1.     ; maximum exposure time (hr)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; set astrophys parameters ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  MzV  = 23.0 ; zodiacal light surface brightness (mag/arcsec**2)
  MezV = 22.0 ; exozodiacal light surface brightness (mag/arcsec**2)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;   set wavelength grid    ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF KEYWORD_SET(compute_lam) THEN BEGIN
    lam  = lammin ;in [um]
    Nlam = 1
    WHILE lam LT lammax DO BEGIN
      lam  = lam + lam/res
      Nlam = Nlam +1
    ENDWHILE
    lam    = FLTARR(Nlam)
    lam[0] = lammin
    FOR j=1,Nlam-1 DO BEGIN
      lam[j] = lam[j-1] + lam[j-1]/res
    ENDFOR
  ENDIF
  Nlam = N_ELEMENTS(lam)
  IF KEYWORD_SET(compute_dlam) THEN BEGIN
    dlam = FLTARR(Nlam) ;grid widths (um)
    FOR j=1,Nlam-2 DO BEGIN
      dlam[j] = 0.5*(lam[j+1]+lam[j]) - 0.5*(lam[j-1]+lam[j])
    ENDFOR
    ;widths at edges are same as neighbor
    dlam[0] = dlam[1]
    dlam[Nlam-1] = dlam[Nlam-2]
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;  set quantum efficiency  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  q = FLTARR(Nlam)
  FOR j=0,Nlam-1 DO BEGIN
    IF (lam[j] LE 0.7) THEN BEGIN
      q[j] = 0.9
    ENDIF ELSE BEGIN
      q[j] = 0.9*(1.0 - (lam[j]-0.7)/(1.0-0.7))
    ENDELSE
    IF q[j] LT 0 THEN q[j] = 0.
  ENDFOR

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;  set dark current and read noise ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  De = FLTARR(Nlam)
  De[*] = 1.d-4
  Re = FLTARR(Nlam)
  Re[*] = 0.1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;  angular size of lenslet ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  theta = lammin/1.d6/diam/2.*(180/!DPI*3600.) ;diameter (assumes sampled at ~lambda/2D) (arcsec)
  IF KEYWORD_SET(nir) THEN BEGIN
    theta = DBLARR(Nlam)
    iVIS  = WHERE(lam LE 1.0)
    iNIR  = WHERE(lam GT 1.0)
    theta[iVIS] = lammin/1.d6/diam/2.*(180/!DPI*3600.)

    IF (iNIR NE -1) THEN BEGIN
      theta[iNIR] = 1.0/1.d6/diam/2.*(180/!DPI*3600.)
      q[iNIR]  = 0.90
      Re[iNIR] = 2.

      ; set dark current based on NIR detector properties
      IF ( lammax LE 2.0 ) THEN De[iNIR] = 1.e-3*10.^( (Tdet-120.)*7./100. )
      IF ( lammax GT 2.0 AND lammax LE 4.0 ) THEN De[iNIR] = 1.e-3*10.^( (Tdet-80.)*9./140. )
      IF ( lammax GT 4.0 AND lammax LE 7.0 ) THEN De[iNIR] = 1.e-3*10.^( (Tdet-40.)*11./140. )
      IF ( lammax GT 7.0 ) THEN De[iNIR] = 1.e-3*10.^( (Tdet-30.)*11./70. )

      ; don't let dark current fall below a threshold
      iDe = WHERE(De[iNIR] LT 1.d-3)
      IF (iDe[0] NE -1) THEN De[iDe] = 1.d-3
    ENDIF

  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;      set throughput      ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  T    = DBLARR(Nlam)
  T[*] = Tput
  sep  = r/d*SIN(alpha*!DPI/180.)*!DPI/180./3600.d ; separation in radians
  iIWA = WHERE( IWA*lam/diam/1.d6 GT sep )
  IF (iIWA[0] NE -1) THEN T[iIWA] = 0. ;points inside IWA have no throughput
  IF KEYWORD_SET(fix_owa) THEN BEGIN
    IF ( OWA*lammin/diam/1.d6 LT sep ) THEN T[*] = 0. ;points outside OWA have no throughput
  ENDIF ELSE BEGIN
    iOWA = WHERE( OWA*lam/diam/1.d6 LT sep  )
    IF (iOWA[0] NE -1) THEN T[iOWA] = 0. ;points outside OWA have no throughput
  ENDELSE

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; degrade albedo spectrum  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF KEYWORD_SET(compute_lam)  THEN A  = DEGRADE_SPEC(Ahr,lamhr,lam,DLAM=dlam)
  IF ~KEYWORD_SET(compute_lam) THEN A  = Ahr
  IF KEYWORD_SET(compute_lam)  THEN Fs = DEGRADE_SPEC(solhr,lamhr,lam,DLAM=dlam)
  IF ~KEYWORD_SET(compute_lam) THEN Fs = solhr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;      compute fluxes      ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;Fs = Fstar(lam, Teff, Rs, r, /AU) ;stellar flux on planet
  Fp = Fplan(A, Phi, Fs, Rp, d)     ;planet flux at telescope
  Cratio = FpFs(A, Phi, Rp, r)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;    compute count rates   ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  cp     =  cplan(q, fpa, T, lam, dlam, Fp, diam)
  cz     =  czodi(q, X, T, lam, dlam, diam, MzV)
  cez    = cezodi(q, X, T, lam, dlam, diam, r, Fstar(lam,Teff,Rs,1.,/AU), Nez, MezV)
  csp    = cspeck(q, T, C, lam, dlam, Fstar(lam,Teff,Rs,d), diam)
  cD     =  cdark(De, X, lam, diam, theta, DNhpix)
  cR     =  cread(Re, X, lam, diam, theta, DNHpix, Dtmax)
  cth    =  ctherm(q, X, lam, dlam, diam, Tsys, emis)
  ctot   = cp + cz + cez + csp + cD + cR + cth
  cnoise = cp + 2*(cz + cez + csp + cD + cR + cth) ; assumes background subtraction

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;  exposure time to SNR=1  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  DtSNR = DBLARR(Nlam)
  DtSNR[*] = 0.d
  i     = WHERE(cp GT 0)
  IF (i[0] NE -1) THEN DtSNR[i] = cnoise[i]/cp[i]^2./3600. ; (hr)

END
