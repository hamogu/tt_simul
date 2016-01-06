;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_RANKINE_HUGENIOT_TURHO
;
; PURPOSE:
;	This procedure wrapps the Rankine-Hugoniot jump conditions
;	for a different set of independent variables
;	without any approximation for weak or strong shocks
;
; CATEGORY:
;	Hydrodynamic
;
; CALLING SEQUENCE:
;	TT_RANKINE_HUGONIOT_TURHO T, U, Rho
;
; INPUTS:
;	T:	Initial temperature
;
;	U:	Initial speed of gas
;
;	Rho:	Initial density
;
; OUTPUTS:
;	The output hs the same format as the input. The results after the shock
;	are stored in the variables in the same units as the input
;
; CALLS:
;	tt_preassure, TT_RANKINE_HUGONIOT, tt_temperature
; RESTRICTIONS:
;	Only for monoatomic gases (gamma=5/3)
;
; TEST STATUS:
;	not tested
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther 21.01.2005
;-


pro TT_RANKINE_HUGONIOT_TURHO, $
	t, $		;temperature
	u, $		;speed of gas
	rho		;density
p=tt_pressure(t,rho)
TT_RANKINE_HUGONIOT, p,u,rho
t=tt_temperature(p,rho)
end
