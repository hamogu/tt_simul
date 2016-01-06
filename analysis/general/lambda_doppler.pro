;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	lambda_doppler
;
; PURPOSE:
;	calculated the doppler shift
;
; CATEGORY:
;	modeling
;
; CALLING SEQUENCE:
;	result=lambda_doppler(lambda_rest,v)
;
; INPUTS:
;	lambda_rest:	rest wavelength 
;	v:		velocity in cm/s
;
; OUTPUTS:
;	result:		relativitically doppler shifted lambda
;
; TEST STATUS:
;	 tested
;
; EXAMPLE:
;	result=lambda_doppler(1032.,300e5)
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 23.11.2005
;-
function lambda_doppler, lambda_rest,v
tt_units
common tt_constants
gamma=(1-(v/c_light)^2)^(-.5)
return, lambda_rest*gamma*(1+v/c_light)
end
