;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_GET_MEAN_ATOMIC_WEIGHT
;
; PURPOSE:
;	This function calculates the mean atomic weight
;
; CATEGORY:
;	CHIANTI-Extension
;
; CALLING SEQUENCE:
;	Result=tt_get_mean_atomic_weight
;
; INPUTS:
;	The abundance is taken from the abundance file specified in !abund
;
; OUTPUTS:
;	The result is given in g
;
; COMMON BLOCKS:
;	Constants:	Needed to get Boltzmanns constant
;	elements:	abundance and ioneq
;
; RESTRICTIONS:
;	assumes usual isotopic abundance
;	see common block constants for further restrictions
;
; CALLS:
;	read_abund
;
; TEST STATUS:
;	run and result seems plausible
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther 20.01.2005
;-


function TT_GET_MEAN_ATOMIC_WEIGHT
common tt_constants
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
if (n_elements(abund) eq 0) || (atomic_mass_unit eq 0) then printf, -2, "common blocks abund or tt_units not initialiezed!"
return, total(abund*mass_number*atomic_mass_unit)/total(abund)
end
