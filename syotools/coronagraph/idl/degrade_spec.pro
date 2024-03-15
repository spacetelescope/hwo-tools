;;;;;;;;;;;;;;;;;;;;   DEGRADE_SPEC   ;;;;;;;;;;;;;;;;;;;;
;    Purpose: Routine to take a spectrum, specified at   ;
;      high spectral resolution, and degrade it to a     ;
;      lower resolution.                                 ;
;    Author: Tyler D. Robinson                           ;
;            robinson@astro.washington.edu               ;
;                                                        ;
;    Inputs:                                             ;
;            specHR - high resolution spectrum, in units ;
;                     of energy/time/area/wavelength     ;
;                     (/solid angle)                     ;
;            lamHR  - high resolution wavelength grid,   ;
;                     units of length must match those   ;
;                     of input spectrum                  ;
;            lamLR  - low resolution wavelength grid,    ;
;                     for output spectrum, same units as ;
;                     input grid                         ;
;                                                        ;
;    Options:                                            ;
;              dlam - if set, contains widths of         ;
;                     spectral elements, assumed to be   ;
;                     centered at lam                    ;
;    Output:                                             ;
;            specLR - low resolution spectrum            ;
;                                                        ;
;    Use:                                                ;
;            specLR =  DEGRADE_SPEC(specHR,lamHR,lamLR)  ;
;                                                        ;
;;;;;;;;;;;;;;;;;;;;   DEGRADE_SPEC   ;;;;;;;;;;;;;;;;;;;;
;
FUNCTION DEGRADE_SPEC, specHR, lamHR, lamLR, DLAM=dlam

  ;;store input variables, so as to not overwrite
  ;;the input quantities
  lamHI  = lamHR
  spec   = specHR
  lamLO  = lamLR
  IF (N_ELEMENTS(dlam) NE 0) THEN dlamLO = dlam
  
  ;;reverse ordering if wavelength vector is
  ;;decreasing with increasing index
  IF lamHR[0] GT lamHR[1] THEN BEGIN
    lamHI = REVERSE(lamHR)
    spec  = REVERSE(specHR)
  ENDIF
  IF lamLR[0] GT lamLR[1] THEN BEGIN
    lamLO  = REVERSE(lamLR)
    IF (N_ELEMENTS(dlam) NE 0) THEN dlamLO = REVERSE(dlam)
  ENDIF

  ;;number of gridpoints in output
  Nspec  = N_ELEMENTS(lamLO)
  ;;declare output spectrum
  specLO = FLTARR(Nspec)
  
  ;;loop over spectral elements
  FOR i=0L, Nspec-1 DO BEGIN
  
    ;;if dlam not provided, we must determine 
    ;;bin widths
    IF (N_ELEMENTS(dlam) EQ 0) THEN BEGIN
      ;;define short and long wavelength edges
      ;;of the low resolution grid, special cases
      ;;at the edge of the grid
      IF (i GT 0 AND i LT Nspec-1) THEN BEGIN
        lamS = 0.5*(lamLO[i-1] + lamLO[i])
        lamL = 0.5*(lamLO[i+1] + lamLO[i])
      ENDIF ELSE IF (i EQ 0) THEN BEGIN
        lamS = lamLO[i] - 0.5*(lamLO[i+1] - lamLO[i])
        lamL = 0.5*(lamLO[i+1] + lamLO[i])
      ENDIF ELSE BEGIN
        lamS = 0.5*(lamLO[i-1] + lamLO[i])
        lamL = lamLO[i] + 0.5*(lamLO[i-1] + lamLO[i])
      ENDELSE
    ENDIF ELSE BEGIN
      lamS = lamLO[i] - 0.5*dlamLO[i]
      lamL = lamLO[i] + 0.5*dlamLO[i]
    ENDELSE

    ;;elements of high resolution grid within 
    ;;the low resolution element
    is = WHERE(lamHI GE lamS AND lamHI LE lamL)
    
    ;;if there aren't any high resolution elements
    ;;within the low resolution element, then error
    IF is[0] EQ -1 THEN PRINT, "Error in DEGRADE_SPECTRUM: wavelength grids do not sync"
    
    ;;if three (3) or less elements of spectrum within 
    ;;the low resolution gridpoint, then do an interpolation, 
    ;;otherwise, integrate the high resolution spectrum over the
    ;;low resolution element, distributing the integrated energy
    ;;into the low resolution element
    IF N_ELEMENTS(is) LE 3 THEN BEGIN
      specs = INTERPOL(spec[is], lamHI[is],lamLO[i])
    ENDIF ELSE BEGIN
      speci = [INTERPOL(spec[is],lamHI[is],lamS),spec[is],INTERPOL(spec[is],lamHI[is],lamL)]
      lami  = [lamS,lamHI[is],lamL]
      specs = INTEGRAL(lami,speci)/(lamL - lamS)
    ENDELSE
  
    ;;insert result into output variable
    specLO[i] = specs
  ENDFOR

  RETURN, specLO

END