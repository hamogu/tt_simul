;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_VORONOV
;
; PURPOSE:
;	This function calculates the dielectronic recombination rate
;
; CATEGORY:
;	DATA FIT
;
; CALLING SEQUENCE:
;	TT_DIELEC_RECOMB, temp,recomb
;
; INPUTS:
;	temp:	electron temperature in K 
;
; OUTPUTS:
;	recomb:		array of recombination rate cm^3/s
;	
; COMMON BLOCKS:
;	elements:	abundance and ioneq
;
; RESTRICTIONS:
;	coeffeciants from H to Ni
;
; NOTES:
;	heavily adapted from proton_dens from the CHIANTI libary
;
; BIBLOGRAPHY:
;	Mazzotta P., Mazzitelli G., Colafrancesco S., Vittorio N.
;        <Astron. Astrophys. Suppl. Ser. 133, 403 (1998)>
;
; TEST STATUS:
;	results seem plausible
;
; MODIFICATION HISTORY:
;	Written by:	Moritz Günther 18.04.2005
;-
pro tt_dielec_recomb,$
		t_e,$		;electron temperature in K
		recomb		;array containing the dielectronic recombination rate in cm^3/s
		
tt_read_array, concat_dir(!tt_data, 'dielecrecomb.dat'), coeff,ref

temp=t_e*8.617385e-05		;from Kelvin in eV
for i=0,3 do recomb[0:27,0:27]+=temp^(-1.5)*coeff[*,*,i]*exp(-coeff[*,*,i+4]/temp)


return
end
