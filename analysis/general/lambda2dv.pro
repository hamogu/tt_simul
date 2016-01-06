;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	lambda2dv
;
; PURPOSE:
;	transforms lambdas in doppler velocities - useful for plot labelling
;
; CATEGORY:
;	modeling
;
; CALLING SEQUENCE:
;	result=lambda2dv(lambda_0,wave)
;
; INPUTS:
;	lambda_0:	centre wavelength 
;	wave:		array of wavelength in Angstroem
;
; OUTPUTS:
;	result:		dv scale in cm/s
;
; TEST STATUS:
;	 tested
;
; EXAMPLE:
;	result=lamda2dv(1032.,findgen(200)/10.+1022.)
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 23.11.2005
;	23.03.2007	Corrected header - M. Guenther
;	23.07.2009	changed signature - blue-shifted lambda is now negative
;-
function lambda2dv, lambda_0,wave
tt_units
common tt_constants
if n_params() ne 2 then message, 'IDL> result=lambda2dv(lambda_0,wave)'
return, double((wave-lambda_0))/double(lambda_0)*c_light
end


