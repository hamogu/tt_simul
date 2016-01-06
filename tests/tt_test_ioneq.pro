;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_TEST_IONEQ
;
; PURPOSE:
;	This procedure computes the ionisation equilibrium from recombination and ionisation rates.
;	Use this for camparison to the literature or the the result of the relaxation test in
;	tt_test_ioneq_rate
;
; CATEGORY:
;	VERIFICATIOM
;
; CALLING SEQUENCE:
;	TT_test_ioneq,result
;      
; OUTPUTS:
;	result:		output an ioneq array as result of simulation in standart CHIANTI form
;			for temperatures from logT=4.0 to 8.0 in 0.1 steps
;
; KEYWORDS:
;	quiet:		suppresses imformational output
;
; RESTRICTIONS:
;	results are only valid as log as the routines called for the rates, e.g. only for iz<29 for cfit.f
;                        
; EXAMPLE:
;	TT_test_ioneq,result
;
; CALLING:
;	tt_dielec_recomb,tt_rad_recomb,tt_arnaudraymond
;
; EXPLANATION:
;	the rate equation can be written in the form Ax=b
;	where A is a tridiagonal matrix of rates, x the ionisation state and b the change dx/dt
;	in equilibrium b=0
;	One of the equations is redundant because the recombination of one ion equals alway the ionisation of another.
;	So the first row is replaced with a normalisation condition: x1+x2+x3+...=1
;
; TEST STATUS:  
;	tested and used
;
; MODIFICATION HISTORY:
;	Written by:	Moritz Günther 28.4.2005
;	5.7.05 Moritz Günther 		added /quiet keyword
;-

pro tt_test_ioneq, result,quiet=quiet

result=make_array(41,30,31,val=0.)

for temp=40,80 do begin
  if ~keyword_set(quiet) then print,'Calculating logT=',float(temp)/10.
  recomb=make_array(31,31,val=0.)
  colion=make_array(31,31,val=0.)
  t_e=10.^(float(temp)/10.)
  ;---calculation recombination---
  tt_dielec_recomb, t_e, recomb
  tt_rad_recomb,t_e,recomb
  ;--- calculation ionization ---
  tt_arnaudraymond, t_e,colion

  for iz=0, (29) do begin		;loop over all elements
    dN_dt=make_array(iz+2,val=0.)	;no change in populations in equilibrium
    ionrate=make_array(iz+2,iz+2,val=0.)
    recombrate=make_array(iz+2,iz+2,val=0.)
    ;fill diagonal elements
    for ion=0,iz+1 do begin		;element Z has Z+1 ionisation stages
      ionrate[ion,ion]=-colion[iz,ion]
      recombrate[ion,ion]=-recomb[iz,ion]
    endfor
    ionrate=ionrate+recombrate-shift(ionrate,0,1)-shift(recombrate,0,-1)	;ionrate is now the total rate
    ;Normalisation: x1+x2+x3+...=1
    ionrate[*,0]=1.
    dn_dt[0]=1.
    ; Decompose ionrate: 
    LUDC, ionrate, INDEX ; Compute the solution using back substitution: 
    result[temp-40,iz,0:iz+1]= LUSOL(ionrate, INDEX, dn_dt)
  endfor  
endfor
end
    