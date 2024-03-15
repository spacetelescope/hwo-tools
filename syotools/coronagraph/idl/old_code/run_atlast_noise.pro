PRO RUN_ATLAST_NOISE, planet=planet, startype=startype, lammin=lammin, lammax=lammax, res=res, Tput=Tput, diam=diam, Tsys=Tsys

 ;planet choices: earth, venus, archean,
                  ;earlymars, hazyarchean, earlyvenus,
                  ;jupiter, saturn, uranus, neptune, mars,
                  ;fstarozone

 ;star choices (for now): sun, f2v
  
  if not keyword_set(planet) then whichplanet = 'earth'
  if keyword_set(planet) then whichplanet = planet

  ;;; planet parameters ;;
  
  alpha = 90.     ;phase angle at quadrature
  Phi   = 1/!DPI  ;phase function at quadrature

  if whichplanet eq 'earth' then begin
     fn = 'earth_quadrature_radiance_refl.dat'
     READCOL, 'planets/'+fn, lamhr, radhr, solhr, reflhr
     Ahr   = !PI*(!PI*radhr/solhr) ;hi-resolution reflectivity
    Rp    = 1.0     ;Earth radii
    r     = 1.0     ;semi-major axis (AU) 
  endif

  if whichplanet eq 'venus' then begin
     fn = 'Venus_geo_albedo.txt'
     READCOL, 'planets/'+fn, lamhr, Ahr
    Rp    = 0.95     ;Earth radii
    r     = 0.72     ;semi-major axis (AU)
  endif

  if whichplanet eq 'archean' then begin
     fn = 'ArcheanEarth_geo_albedo.txt'
     READCOL, 'planets/'+fn, lamhr, Ahr
     Rp    = 1.0     ;Earth radii
     r     = 1.0     ;semi-major axis (AU)
  endif

  if whichplanet eq 'earlymars' then begin
     fn = 'EarlyMars_geo_albedo.txt'
     READCOL, 'planets/'+fn, lamhr, Ahr
     Rp    = 0.53     ;Earth radii
     r     = 1.52     ;semi-major axis (AU)
  endif

  if whichplanet eq 'hazyarchean' then begin
     fn = 'Hazy_ArcheanEarth_geo_albedo.txt'
     READCOL, 'planets/'+fn, lamhr, Ahr
     Rp    = 1.0     ;Earth radii
     r     = 1.0     ;semi-major axis (AU)
  endif

  if whichplanet eq 'earlyvenus' then begin
     fn = 'EarlyVenus_geo_albedo.txt'
     READCOL, 'planets/'+fn, lamhr, Ahr
    Rp    = 0.95     ;Earth radii
    r     = 0.72     ;semi-major axis (AU)
  endif

  if whichplanet eq 'jupiter' then begin
     fn = 'Jupiter_geo_albedo.txt'
     READCOL, 'planets/'+fn, lamhr, Ahr
    Rp    = 10.86     ;Earth radii
    r     = 5.20     ;semi-major axis (AU)
  endif

  if whichplanet eq 'saturn' then begin
     fn = 'Saturn_geo_albedo.txt'
     READCOL, 'planets/'+fn,lamhr, Ahr
    Rp    = 9.00     ;Earth radii
    r     = 9.54     ;semi-major axis (AU)
  endif

  if whichplanet eq 'uranus' then begin
     fn = 'Uranus_geo_albedo.txt'
     READCOL, 'planets/'+fn, lamhr, Ahr
    Rp    = 3.97     ;Earth radii
    r     = 19.19     ;semi-major axis (AU)
 endif


  if whichplanet eq 'warmuranus' then begin
     fn = 'Uranus_geo_albedo.txt'
     READCOL, 'planets/'+fn, lamhr, Ahr
    Rp    = 3.97     ;Earth radii
    r     = 5.20     ;semi-major axis (AU)
 endif

  if whichplanet eq 'warmneptune' then begin
     fn = 'Neptune_geo_albedo.txt'
     READCOL, 'planets/'+fn, lamhr, Ahr
    Rp    = 3.97     ;Earth radii
    r     = 5.20     ;semi-major axis (AU)
 endif

  if whichplanet eq 'neptune' then begin
     fn = 'Neptune_geo_albedo.txt'
     READCOL, 'planets/'+fn, lamhr, Ahr
    Rp    = 3.85     ;Earth radii
    r     = 30.07     ;semi-major axis (AU)
 endif

  if whichplanet eq 'mars' then begin
     fn = 'Mars_geo_albedo.txt'
     READCOL, 'planets/'+fn, lamhr, Ahr
    Rp    = 0.53     ;Earth radii
    r     = 1.52     ;semi-major axis (AU)
 endif

  if whichplanet eq 'fstarozone' then begin
     fn = 'fstarcloudy.sav'
     ;fn = 'F2V_5.e-1fCO2_1.e6H2Volc_1.e10BIF.out_toa.rad'
     ;READCOL, 'planets/'+fn, lamhr, wnhr, sunhr, fluxhr
     restore, 'planets/'+fn
     lamhr = reverse(wl_)
     ahr = reverse(cloudyalb_)
     ;Ahr = (2/3.) * fluxhr/(sunhr/2.) ;geometric albedo & phase corrections
       Rp    = 1.0     ;Earth radii
       r     = 1.72    ;semi-major axis (AU)
  endif



  ;;; star parameters ;;;
  if not keyword_set(startype) then begin
  Teff  = 5780.   ;Sun-like Teff (K)
  Rs    = 1.      ;star radius in solar radii
endif
  if keyword_set(startype) then begin
  if  startype eq 'sun' then begin
  Teff  = 5780.   ;Sun-like Teff (K)
  Rs    = 1.      ;star radius in solar radii
endif
  if  startype eq 'f2v' then begin
  Teff  = 7050.   ;Sun-like Teff (K)
  Rs    = 1.3      ;star radius in solar radii
endif
endif

  
  ;;; system parameters ;;;
  d    = 10.     ;distance to system (pc)
  Nez  = 1.      ;number of exo-zodis
 
  
  ;;; telescope parameters ;;;
  if not keyword_set(diam) then diam   = 8.     ;diameter (m)
  C      = 1.d-10  ;raw contrast
  IF not keyword_set(lammin) then lammin = 0.3     ;minimum wavelength (um)
  IF not keyword_set(lammax) then  lammax = 2.0     ;maximum wavelength (um)
  IF not keyword_set(Res) then Res    = 200.     ;spectral resolution (lambda/Dlambda)
  IF not keyword_set(Tput) then Tput   = 0.2     ;system throughput
  IWA    = 0.5      ;inner working angle (lambda/D)
  OWA    = 30000.     ;outer working angle (lambda/D)
  if not keyword_set(Tsys) then Tsys   = 180.    ;telescope/system temperature (K)
  emis   = 0.9     ;telescope/system emissivity
  wantsnr = 10.     ;desired SNR
  
  ;Shawn: "I don't like noise.  It makes me sad."

  ATLAST_NOISE, Ahr, lamhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez, $
                lammin, lammax, Res, diam, Tput, C, IWA, OWA, Tsys, emis, $
                lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, $
                cth, DtSNR, wantsnr, whichplanet, /COMPUTE_LAM


;integration 
;R = 70
;3 curves: throughput 10%, 
;8 m 20%
;4 m 20%
;4 m 50% 


END
