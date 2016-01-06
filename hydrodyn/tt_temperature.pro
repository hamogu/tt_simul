;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_TEMPERATURE
;
; PURPOSE:
;	This function calculates the temperatur for an ideal gas
;	give the temperature and the density
;
; CATEGORY:
;	Hydrodynamic
;
; CALLING SEQUENCE:
;	Result=TT_TEMPERATURE(p,Rho)
;
; INPUTS:
;	p:	preassure in dyn/cm^2
;
;	Rho:	Density in g/cm^(-3)
;
; OUTPUTS:
;	The result is given in Kelvin
;
; COMMON BLOCKS:
;	TT_Constants:	Needed to get Boltzmanns constant
;	TT_mean_atomic_weight: 	Well, obviously here for mean_atomic_weight
;
; FUNCTION:
;	temp=p*mean_atomic_weight/(k_boltz*rho)
;
; RESTRICTIONS:
;	Only for ideal gas (PV=NkT)
;
; TEST STATUS:
; 	tested with real data
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther 20.01.2005
;-


function TT_TEMPERATURE, $
	p, $		;pressure
	rho		;density
common TT_constants
common TT_mean_atomic_weight
return, p*mean_atomic_weight/(k_boltz*rho)
end
