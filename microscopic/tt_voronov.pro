;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_VORONOV
;
; PURPOSE:
;	This function calculates the ionisation rate
;
; CATEGORY:
;	DATA FIT
;
; CALLING SEQUENCE:
;	TT_VORONOV, temp,colion
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
; NOTES:
;	heavily adapted from proton_dens from the CHIANTI libary
;
; BIBLOGRAPHY:
;	G. S. Voronov, 1997, ADNDT, 65, 1
;	translated to IDL from cfit(iz,in,t,c) Version 2, March 24, 1997
;	at http://www.pa.uky.edu/~verner/fortran.html
;
; TEST STATUS:
;	partly tested
;
; MODIFICATION HISTORY:
;	Written by:	Moritz Günther 31.01.2005
;-
pro tt_voronov,$
		t_e,$		;electron temperature in K
		colion		;array containing the collisinal ionisation rate in cm^3/s
		
tt_read_array, concat_dir(!tt_data, 'voronov.dat'), coeff,ref
coeff=coeff[1:*,1:*,*]		;cut out the 0 th line and row
u=coeff[*,*,0]/(t_e*8.617385e-05)
; formula is invalid for cases with u gt 80
u=u-(u gt 80.)*u
s=size(u)
matrix1=make_array(s[1:2],value=1.)

colion=coeff[*,*,2]*(matrix1+coeff[*,*,1]*sqrt(u))/(coeff[*,*,3]+u)*u^coeff[*,*,4]*exp(-u)
;set 0 where non applicable and therefor the above calculation may
;result in a "-NaN"
colion[where(u eq 0)]=u[where( u eq 0)]


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
;check this statement!
for j=0,n_elements(colion[*,0])-1 do colion[j,0:j]=reverse(reform(colion[j,0:j]),1)
return
end
