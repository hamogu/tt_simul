;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_GET_XE
;
; PURPOSE:
;	This function calculates the electron to (atom+ion) ratio
;
; CATEGORY:
;	CHIANTI-Extension
;
; CALLING SEQUENCE:
;	Result=TT_GET_XE(temp)
;
; INPUTS:
;	temp:	Temperature in K 
;
; OUTPUTS:
;	ratio electron/(atom+ion)=xe
;
; COMMON BLOCKS:
;	elements:	abundance and ioneq
;
; NOTES:
;	heavily adapted from proton_dens from the CHIANTI libary
;	works only if ioneq is cut to a single temperature beforehand
;
;
; TEST STATUS:
;	tested
;
; MODIFICATION HISTORY:
; 	Adapted from:	proton_dens V.4
;	Written by:	Moritz Günther 27.01.2005
;-


function tt_get_xe,$
	temp_in		;Temperature in K

common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

siz=size(ioneq)
n_ion=siz[3]

i_a=where(abund NE 0.,nia)
nelec=0d0
 
FOR i=0,nia-1 DO BEGIN
  nelec=nelec+abund[i_a[i]]*(transpose(reform(ioneq[0,i_a[i],*])) # findgen(n_ion))
ENDFOR
   
return,nelec[0]/total(abund)

END
