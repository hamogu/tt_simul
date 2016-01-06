;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	PREPARE_ANALYSIS_HYD
;
; PURPOSE:
;	This procedure calculates the basic diagnostcs of a plasma from a simulated hydrodynamical structure.
;
; CATEGORY:
;	analysis
;
; CALLING SEQUENCE:
;	prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens
;
; INPUTS:
;	hydrodyn:	hydrodynmic shock structure as output by tt_simul
;

;
; OUTPUTSS:
;	delta_x:	a vector with the stepsize between different point in hydrodyn
;	volem:		a vector containing the volume emission meassure of every step [cm^{-3}]
;	elecdens:	a vector containing the electron density in every step [cm^{-3}]
;
; TEST STATUS:
;	checked
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 22.08.2005
;	transformed in a procedue for hydrodynamics:	Moritz Guenther 3.11.2005
;-

pro prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
n_steps=n_elements(hydrodyn[0,*])
delta_x=-hydrodyn[0,*]+[[reform(max(hydrodyn[0,*]),1,1)],[hydrodyn[0,0:n_steps-2]]]
delta_te=-hydrodyn[5,*]+[[reform(max(hydrodyn[5,*]),1,1)],[hydrodyn[5,0:n_steps-2]]]
volem=delta_x*hydrodyn[2,*]^2*hydrodyn[6,*]
elecdens=hydrodyn[2,*]*hydrodyn[6,*]
end