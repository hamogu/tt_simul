;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_PRESSURE
;
; PURPOSE:
;	This function calculates the preassure for an ideal gas
;	give the temperature and the density
;
; CATEGORY:
;	Hydrodynamic
;
; CALLING SEQUENCE:
;	Result=TT_PRESSURE(T,Rho)
;
; INPUTS:
;	T:	Temperature in Kelvin
;
;	Rho:	Density in g7cm^(-3)
;
; OUTPUTS:
;	The result is given in dyn/cm^2
;
; COMMON BLOCKS:
;	tt_Constants:	Needed to get Boltzmanns constant
;	tt_mean_atomic_weight: 	Well, obviously here for mean_atomic_weight
;
; FUNCTION:
;	p=temp*K_boltz*rho/mean_atomic_weight
;
; RESTRICTIONS:
;	Only for ideal gas (PV=NkT)
;
; TEST STATUS:
;	tested with real data
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther 20.01.2005
;-


function TT_PRESSURE, $
	t, $		;temperature
	rho		;density
common tt_constants
common tt_mean_atomic_weight
return, t*k_boltz*rho/mean_atomic_weight
end
