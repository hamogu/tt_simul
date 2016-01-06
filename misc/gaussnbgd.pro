;+
; NAME:
;	GAUSSNBGD
;
; PURPOSE:
;	Returns values (and partial derivatives) of a sum of gaussian with a const background.
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
;	params:		param[0] is const background, the rest is for gaussn
;
; OPTIONAL INPUTS:
;	partialderiv:	a named variable containing the partial derivatives on output
;
; OUTPUT:
;	Value of function at input x (possibly array)
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 26.01.2007
;-


pro gaussnbgd,x,params,f,partialderiv
case (n_elements(params) mod 3) of
  1: begin
     gaussn,x,params[1:*],f,partialderiv
     f+=params[0]
     if n_elements(partialderiv) ne 0 then partialderiv=[[make_array(n_elements(x),val=1.)],[partialderiv]]
     end
  else: message,'Not implementd yet!'
endcase
end