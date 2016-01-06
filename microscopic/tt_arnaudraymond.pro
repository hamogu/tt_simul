;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_ArnaudRaymond
;
; PURPOSE:
;	This function calculates the ionisation rate
;
; CATEGORY:
;	DATA FIT
;
; CALLING SEQUENCE:
;	TT_ArnaudRaymond, temp,colion
;
; INPUTS:
;	temp:	electron temperature in K 
;
; OUTPUTS:
;	colion:		array of colionination rate cm^3/s
;	
; COMMON BLOCKS:
;	elements:	abundance and ioneq
;
;CALLING:
;	colfit.c: f2c translated into C and modified to match the IDL protable calling convention
;	see documentation of colfit.c
;
; BIBLOGRAPHY:
;	M. Arnaud and J. Raymond, 1992, ApJ 398, 394 for Fe 
;	M. Arnaud and R. Rothenflug, 1985, A&AS 60, 425 for H, He, C, N, O, Ne,
;	Na, Mg, Al, Si, S, Ar, Ca, Ni (with some corrections described in
;	Verner & Yakovlev, 1990, ApSS, 165, 27) 
;	interpolation/extrapolation for other elements.
;	Fortran code 
;	Version 4. January 8, 1997.
;	at http://www.pa.uky.edu/~verner/fortran.html
;
; TEST STATUS:
;	checked for plausibility and compared to tt_voronov
;
; MODIFICATION HISTORY:
;	Written by:	Moritz Günther 02.03.2005
;-
pro tt_arnaudraymond,$
		t_e,$		;electron temperature in K
		colion		;array containing the collisinal ionisation rate in cm^3/s
		
replicate_inplace, colion, 0.	;faster than colion[*]=0.

path=concat_dir(!tt_data,'colfit.so')
FOR iz=1,n_elements(colion[*,0]) DO BEGIN		;n_elements(recomb[*,0]) = number of chemical elements
    FOR ion=1,iz DO BEGIN
        colion[iz-1,iz-ion]=call_external(path,"colfit",iz,ion,0,float(t_e), /F_VALUE)
    ENDFOR
ENDFOR
;Problem is the ordering of the input data:
;(HI	HeII	LiIII	...)
;(0	HeI	LiII	...)
;(0	0	LiI	...)
;(...			   )
;need to change that to standard in ioneq
;(HI	HeI	LiI	...)
;(0	HeII	LiII	...)
;(0	0	LiIII	...)
;(...			   )
;that is why filling colion[iz-1,iz-ion]


return
end
