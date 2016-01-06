;+
; NAME:
;	GAUSSN
;
; PURPOSE:
;	Returns values (and partial derivatives) of a sum of gaussian.
;	Meant for use with curve_fit
;
; CATEGORY:
;	Mathematics
;
; CALLING SEQUENCE:
;	GAUSSN,x,params,f[,partialderiv]
;
; INPUTS:
;	x:		values for which the function shall be computed, possibly array
;	params:		array with n*3 elements, where n is the number of gaussians with the entries
;			i  : Maximum value of gaussian
;			i+1: Centre
;			i+2: sigma    for i=0 to n-1
;
; OUTPUT:
; Value of function at input x (possibly array)
;
; OPTIONAL OUTPUTS:
;	partialderiv:	a named variable containing the partial derivatives on output
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 25.01.2007
;-


pro gaussn,x,params,f,partialderiv
if (n_elements(params) mod 3) eq 0 then n=n_elements(params)/3 else message,'%GAUSSN: Give three elements in param array per Gaussian!'
f=make_array(n_elements(x),val=0)
if n_params() eq 3 then begin
  partialderiv=make_array(n_elements(f))
  partder=[0,0,0]
  for i=0,n-1 do begin
    f+=gaussian(x,params[i*3:(i*3)+2],partialderiv[i*3:(i*3)+2])
    partialderiv[i*3:(i*3)+2]=partder
  endfor
endif else begin
  for i=0,n-1 do f+=gaussian(x,params[i*3:(i*3)+2])
endelse
end