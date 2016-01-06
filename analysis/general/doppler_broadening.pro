;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	doppler_broadening
;
; PURPOSE:
;	calculates the doppler broadning
;
; CATEGORY:
;	modeling
;
; CALLING SEQUENCE:
;	result=doppler_broadening(lambda_0,temp,iz)
;
; INPUTS:
;	lambda_0:	centre wavelength 
;	temp:		gas temperature in K
;	iz:		charge of ion (to get mass)
;
; OUTPUTS:
;	result:		relativitically doppler shifted lambda
;
; TEST STATUS:
;	 tested
;
; EXAMPLE:
;	result=doppler_broadening(1032.,3e6,6)
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 23.11.2005
;-
function doppler_broadening,lambda_0,temp,iz
tt_units
common tt_constants
return, lambda_0/c_light*(2*k_Boltz*temp/(mass_number(iz-1)*atomic_mass_unit))^.5
end


