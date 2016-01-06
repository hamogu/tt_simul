;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_RANKINE_HUGENIOT
;
; PURPOSE:
;	This procedure evalutes the Rankine-Hugoniot jump conditions
;	without any approximation for weak or strong shocks
;
; CATEGORY:
;	Hydrodynamic
;
; CALLING SEQUENCE:
;	TT_RANKINE_HUGONIOT P, U, Rho
;
; INPUTS:
;	P:	Initial Preassure
;
;	U:	Initial speed of gas
;
;	Rho:	Initial density
;
; OUTPUTS:
;	The output hs the same format as the input. The results after the shock
;	are stored in the variables in the same units as the input
;
; RESTRICTIONS:
;	Only for monoatomic gases (gamma=5/3)
;
; TEST STATUS:
;	run with test values
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther 20.01.2005
;-


pro TT_RANKINE_HUGONIOT, $
	p, $		;pressure
	u, $		;speed of gas
	rho		;density

gamma=5./3.		;monoatomic gas
p_new=(2*rho*u^2-(gamma-1)*p)/(gamma+1)
u_new=(p-p_new)/(rho*u)+u
rho_new=rho*u/u_new
;set return values
p=p_new
rho=rho_new
u=u_new
end
