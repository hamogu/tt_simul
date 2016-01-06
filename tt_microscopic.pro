;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_MICROSCOPIC
;
; PURPOSE:
;	This procedure does the main task of simulating the microscopic physics.
;	It calculates the ionisation/recombination processes and calls
;	CHIANTI routines to find the radiative loss (in ff,fb and bb processes)
;
; CATEGORY:
;	microscopic
;
; CALLING SEQUENCE:
;	tt_microscopic, t_ion,T_e,n,x_e,d_xe,en_loss,relabundHII,out,radtemp=radtemp,rate_test=rate_test
;
; INPUTS:
;	T_ion:		Ion temperature in Kelvin
;	t_e:		Electron temperature in Kelvin
;	n:		Density in 1/cm^(-3)
;	out:		reference of an output object
;
; OPTIONAL INPUTS:
;	radtemp:	radiation of background blackbody, influences rad_loss via different level populations

; OUTPUTS:
;	x_e:		number of free electrons per ion/atom
;	d_xe:		d(x_e)/dx
;	en_loss:	energy lost from the system (due to radiative processes/ionizations/...) in später
;	relativeabundHII: ratio H-II/total number of heavy particles
;
; KEYWORDS:
;	rate_test:	the time consuming rad_loss calculations are not done
;
; COMMON BLOCKS:
;	tt_Constants:	Needed to get Boltzmanns constant
;	tt_units:	defines units symbols (km,m,cm,...)
;	tt_mean_atomic_weight: 	Well, obviously here for mean_atomic_weight
;	tt_steps:	for communication of step size delta_x and surface element A
;	tt_micros_internal: availability of variables throughout the functions of this module
;	tt_elements:	abund and ioneq - from CHIANTI
;
; ADDITIONAL FUNCTIONS
;	This file consists of the main procedure and a supplementary function
;	Function tt_micro_deriv: calculates the derivative for the ODE solver
;
; ACCURACY
;	lsode should solve the ODEs to an absolute and relative error of 1e-07
;	assuming that the hydordynamic variables do not change
;	due to numerical problems of lsode sometimes even negative values occur close
;	to the numerical accuracy. Values seem reliable a low as abour 1e-4.
;
; TEST STATUS:
;	results give correct qualitative behaviour
;	extensive testing with tt_tests and tt_test_time_evolution succsesful
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 15.05.2005
;	
;-

Function tt_micro_deriv,t,y
  common tt_micros_internal, A
  return, reform(A##y)
end 

PRO TT_MICROSCOPIC, $
	t_ion,$			;ion temperature
	t_e,$			;electron temperature
	n_ion,$			; in 1/cm^3
	x_e,$			;Output: ratio electrons/nuclii
	d_xe,$			;Output: change in x_e
	en_loss,$		;Output: energy loss rate in ergs/s
	relabundHII,$		;
	out,$
	rate_test=rate_test,$
	radtemp=radtemp
common tt_units
common tt_constants
common tt_mean_atomic_weight
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
common tt_steps
common tt_micros_internal, ionrate

;später add check input params
recomb=make_array(n_elements(ioneq[0,*,0]),n_elements(ioneq[0,0,*]),val=0.)
colion=make_array(n_elements(ioneq[0,*,0]),n_elements(ioneq[0,0,*]),val=0.)
recomben=make_array(n_elements(ioneq[0,*,0]),n_elements(ioneq[0,0,*]),val=0.)

;use electron temperature in procedures which read data from common block elements
effcoltemp=t_e*(1+t_ion/t_e*m_e/mean_atomic_weight)	;effective collision temperature
ioneq_logt=alog10(effcoltemp)

n_e=x_e[0]*n_ion		;to get a number not an array if x_e is an array of one element

;---calculation recombination---
tt_dielec_recomb, effcoltemp, recomb
tt_rad_recomb,effcoltemp,recomb

;--- calculation ionization ---
tt_arnaudraymond, effcoltemp,colion
n_elem=n_elements(colion[*,0])<n_elements(ioneq[0,*,0])
; size of colion depends on the number of coeff in interpolation
; and could be different from ioneq
colion=colion[0:(n_elem-1),0:(n_elem)]*n_e[0]
recomb=recomb[0:(n_elem-1),0:(n_elem)]*n_e[0]
;colion is in cm^3/s

;calculating energy loss due to ionization
read_ip,!xuvtop+'/ip/chianti.ip',ionpot,ipref
ion_pot=ionpot*cm*c_light*h_planck*e_charge*2.*!pi 	;convert energy from cm^(-1)in erg

;--- calculation of new ionisation state ---

en_loss=0.

;There are different methods to solve the system of ODEs describing the ionisation state
;but the numerical method works best
for iz=0, (n_elem-1) do begin		;loop over all elements
  if (total(ioneq[0,iz,*]) gt 0.) then begin
    ;--- build A matrix for y'=A y---
    ionrate=make_array(iz+2,iz+2,val=0.)
    recombrate=make_array(iz+2,iz+2,val=0.)
    ;fill diagonal elements
    for ion=0,iz+1 do begin		;element Z has Z+1 ionisation stages
      ionrate[ion,ion]=-colion[iz,ion]
      recombrate[ion,ion]=-recomb[iz,ion]
    endfor
    ionrate=ionrate+recombrate-shift(ionrate,0,1)-shift(recombrate,0,-1)	;ionrate is now the rate for the ionisation in y'=A*y
    oldioneq=reform(ioneq[0,iz,0:iz+1])		;ionisation states of element iz
    ;---solve using lsode---
    status=1
    newioneq=lsode(oldioneq,0,delta_t,'tt_micro_deriv',status)
    if (status lt 1) then out->plog,'Numerical problems in integration, status of lsode:'+string(status)
    if (status lt -1) then out->plog,'ERROR in numerical integration'
    ;values of newioneq may get small negative numbers within numerical precission
    ;too high precission wastes time and could be problematic in floating numbers
    if (min(newioneq) lt -1e-7) then begin
      out->log,'For element'+string(iz+1)
      out->plog,'negative value in ioneq:'+string(min(newioneq))
    endif  
    ;the order of magnitude of the values is often correct
    newioneq=abs(newioneq)
    ;---energy losses due to ionization and recombination
    
    ;if change in ioneq of one ionisation state is less then 0.2% of all ions then the error of energy loss
    ;by using just the final value is less then 0.1% which is accepted
    if (abund[iz]*max(abs(newioneq-oldioneq))/2 gt 0.001) then out->plog, 'Large change ioneq - relative inaccuracy in radiative energy loss may be larger than:'+string(abund[iz]*max(newioneq-ioneq[0,iz,0:iz+1])/2)  
    
    en_loss+=total(abund[iz]*newioneq*(ion_pot[iz,0:iz]*colion[iz,0:iz])); +recomben[iz,0:iz+1]*recomb[iz,0:iz+1])
    
    ioneq[0,iz,0:iz+1]=newioneq
    test=check_math(mask=32)	;will reset math error status (32 is floating underflow)
  endif
endfor  

for iz=0,n_elem-1 do begin  
  if (total(ioneq[0,iz,0:n_elem-1]) gt 0) then begin 
    ;ensure ioneq is normalized and cannot run away due to numerical errors
    ;but only for elements with non zero ioneq values
    ;check if one row departs further than a bit from 1
    if (total(ioneq[0,iz,0:n_elem]) lt 0.95 or total(ioneq[0,iz,0:n_elem]) gt 1.05) then out->plog,'Large numerical Error in ioneq calculation - check cause!'
    ioneq[0,iz,0:n_elem]/=total(ioneq[0,iz,0:n_elem])
  endif
endfor  

if not keyword_set(rate_test) then begin
  ; --- calculate radiative energy losses ---
  ;free-free
  ff_rad_loss, t,q_ff_rad_loss, /no_setup
  if (n_elements(t)gt 1) then out->log, "ERROR in tt_microscopic.pro, more than 1 element in ioneq_logt"
  ;bound-bound
  ;no proton rates because T_ion ne T_e, need to change CHIANTI code first
  bb_rad_loss, t,q_bb_rad_loss, /no_setup, /noprot, density=n_e,radtemp=radtemp,rphot=1.
  ;free-bound
  tt_fb_rad_loss_variante2,q_fb_rad_loss,t
  en_loss+=q_ff_rad_loss[0]+q_fb_rad_loss[0]+q_bb_rad_loss[0] 	;rad_loss[0]+ion_loss[0]
endif
if not (where(ioneq lt 0.) eq -1) then out->plog, 'ERROR - Negative Value in ioneq!'
x_e_old=x_e
x_e=tt_get_xe(t_ion)
;---output---
;relative abundance of HII
relabundHII=abund[0]/total(abund)*ioneq[0,0,1]
d_xe=(x_e-x_e_old)/delta_x

;There are always floating-point underflows in this procedure, no reason to worry

end
