;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_ArnaudRaymond
;
; PURPOSE:
;	This function calculates the radiative recombination rate
;
; CATEGORY:
;	DATA FIT
;
; CALLING SEQUENCE:
;	TT_RAD_RECOMB, temp,recomb
;
; INPUTS:
;	temp:	electron temperature in K 
;
; OUTPUTS:
;	recomb:		array of radiative recombination rate cm^3/s
;			added to recomb on input
;	
; COMMON BLOCKS:
;	elements:	abundance and ioneq
;
; CALLING:
;	rrfit.c: f2c translated into C and modified to match the IDL protable calling convention
;	see documentation of rrfit.c
;
; RESTRICTIONS:
;	works from H to Zn
;
; BIBLOGRAPHY:
;	H-like, He-like, Li-like, Na-like - Verner & Ferland, 1996, ApJS, 103, 467
;	Other ions of C, N, O, Ne - Pequignot et al. 1991, A&A, 251, 680,
;	   refitted by Verner & Ferland formula to ensure correct asymptotes
;	Fe XVII-XXIII - Arnaud & Raymond, 1992, ApJ, 398, 394
;	Fe I-XV - refitted by Verner & Ferland formula to ensure correct asymptotes
;	Other ions of Mg, Si, S, Ar, Ca, Fe, Ni - Shull & Van Steenberg, 1982, ApJS, 48, 95
;	Other ions of Na, Al - Landini & Monsignori Fossi, 1990, A&AS, 82, 229
;	Other ions of F, P, Cl, K, Ti, Cr, Mn, Co (excluding Ti I-II, Cr I-IV,
;	Mn I-V, Co I)        - Landini & Monsignori Fossi, 1991, A&AS, 91, 183
;	All other species    - interpolations of the power-law fits
;	Fortran code 
;	Version 4. June 29, 1999.
;	at http://www.pa.uky.edu/~verner/fortran.html
;                        
; TEST STATUS:  
;	no tests
;
; MODIFICATION HISTORY:
;	Written by:	Moritz Günther 18.4.2005
;-
pro tt_rad_recomb,$
		t_e,$		;electron temperature in K
		recomb		;array containing the radiative recombination rate in cm^3/s
		

path=concat_dir(!tt_data,'rrfit.so')
FOR iz=1,(n_elements(recomb[*,0])<30) DO BEGIN		;n_elements(recomb[*,0]) = number of chemical elements
    FOR n_elec=1,iz DO BEGIN
        recomb[iz-1,iz-(n_elec-1)]+=call_external(path,"rrfit",iz,n_elec,float(t_e), /F_VALUE)
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
;that is why filling recomb[iz-1,iz-ion]


return
end
