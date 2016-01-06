;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_TEST_RATE
;
; PURPOSE:
;	test tt_dielec_recomb, tt_rad_recomb and tt_arnaudraymond
;	loops over all temperatures
;
; CATEGORY:
;	VERIFICATIOM
;
; CALLING SEQUENCE:
;	TT_test_ioneq_rate,rec,col
;      
; OUTPUTS:
;	rec:	array of recombination coefficients
;	col:	array of ionision coefficients
;                        
; EXAMPLE:
;	tt_test_rate, rec,col
;
; TEST STATUS:  
;	tested and used
;
; MODIFICATION HISTORY:
;	Written by:	Moritz Günther 21.4.2005
;-
pro TT_test_rate,rec,col 

common tt_units
common tt_steps, delta_x, A, delta_t,x,stepcoeff



;--- Setup variables and calculate starting condition behind shock ---
print, "TT_SIMUL"
print, "Simulation of spectra from the accreting region of T Tauri stars"
print, ""
print, "Moritz Günther"
print, "Hamburger Sternwarte"
print, ""
a=1*cm^2
;initialize output


delta_x=1e+4		;set any non-zero starting value
delta_t=1e0*s

max_number_of_steps=30
rec=make_array(41,31,31)
col=make_array(41,31,31)



for temp=40,80 do begin
  t_e=10.^(float(temp)/10.)
  recomb=make_array(31,31,val=0.)
  colion=make_array(31,31,val=0.)
  ;---calculation recombination---
  tt_dielec_recomb, t_e, recomb
  tt_rad_recomb,t_e,recomb
  ;--- calculation ionization ---
  tt_arnaudraymond, t_e,colion
  rec[temp-40,*,*]=recomb
  col[temp-40,*,*]=colion
endfor

end


