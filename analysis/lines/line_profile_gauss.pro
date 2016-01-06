;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	line_profile_gauss
;
; PURPOSE:
;	this procedure produces a normilised Gauss distribution for line profile calculations
;
; CATEGORY:
;	modeling
;
; CALLING SEQUENCE:
;	result=line_profile_gauss( lambda_0,wave_grid,broadening)
;
; INPUTS:
;	lambda_0:	rest wavelength 
;	wave_grid:	wavelength grid
;	broadening:	width of the gauss in wavelength units
;
; OUTPUTS:
;	result:		an array of the same size as wave_gridwith the guass distribution centered at lamba_0
;
; TEST STATUS:
;	 looks like a great gauss
; EXAMPLE:
;	result=line_profile_gauss(1032.,findgen(200)/10.+1022,1.)
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 23.11.2005
;-
function line_profile_gauss, lambda_0,wave_grid,broadening

common tt_constants

if n_params() lt 3 then begin
   message,"Necessary parameters missing" ,/continue
   return,-1
endif

delta_lambda=lambda_0-wave_grid

if lambda_0 lt min(wave_grid) or lambda_0 gt max(wave_grid) then message,'rest wavelength not included in wavelength interval',/informational
if broadening gt 0 then begin 
  ;result=(!pi)^(-.5)*lambda_0^2/(c_light*broadening)*exp(-(delta_lambda/broadening)^2)
  result=(!pi)^(-.5)/(broadening)*exp(-(delta_lambda/broadening)^2)
endif else begin
  message,'Line braodening with coefficient <=0. Returning no profile', /continue
  result=make_array(n_elements(wave_grid),val=0)
  test=min(abs(wavegrid-lambda_0),min_index)
  result[min_index]=1
endelse
return, result
end

