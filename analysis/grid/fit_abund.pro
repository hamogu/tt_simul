;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	fit_abund
;
; PURPOSE:
;	This procedure calculates abundance ratios. It compares a measured line ratio with a simulated one and the ratio of these ratios is than taken as
;	the element abundance. This is all normalised to the ratios, which are used in the TT_SIMUL files.
;
; CATEGORY:
;	fitting
;
; CALLING SEQUENCE:
;	FIT, ratios, datafile,lev=lev,ratiofile=ratiofile
;
; INPUTS:
;	ratios:		a string array with the ratios to be computed, in a format compatible to calculate_ratios.pro
;	datafile:	path to observational data in a format compatible to read_obs_data.pro
;	index_min:	1-dim index of the best-fit model
;
; OPTIONAL INPUTS:
;	simdata_file:	path to an IDL .sav with with a data structure "data" which contains an output of analyse_grid.pro
;			default is rad6e3newfe.sav in the working directory
;	nh		hydrogen column density in cm^-2 which is used for correction of the observations
;
; KEYWORDS:
;	ratiofile:	if set ratios is not a string array, but the path to a file with such an array written line by line	
;
; DATA:
; 	Uses the data from TT_SIMUL files. The simulation runs with Allen(1975) abundances. This is hardwired in the display code.
;
; TEST STATUS:
;	tested
;
; MODIFICATION HISTORY:
;	Written by:	Moritz Guenther 23.09.2005
;-


pro fit_abund, ratios, datafile,index_min,ratiofile=ratiofile,simdata_file=simdata_file,nh=nh
common tt_units
if n_elements(simdata_file) eq 1 then begin
  restore,simdata_file 
  simdata=data
endif else begin 
  restore,'rad6e3newfe.sav'
  simdata=data
endelse  

read_obs_data,datafile,obs,obs_ref,nh=nh

calculate_ratios, ratios,obs,simdata,obsratios,simratios,ratiofile=ratiofile

n_ratio=n_elements(obsratios[*,0])

abundances=obsratios
abundances[*]=0
;values
abundances[*,0]=obsratios[*,0]/simratios[*,index_min]
;errors
abundances[*,1]=obsratios[*,1]/simratios[*,index_min]
n_ratios=n_elements(abundances[*,0])
abund=total(abundances[*,0]/abundances[*,1]^2.)/total(1./abundances[*,1]^2.)
abund_sigma_gauss=sqrt(1./total(1./(abundances[*,1]^2.)))
abund_sigma_direct=n_ratio le 1 ?0:sqrt(total((abundances[*,0]-abund)^2./abundances[*,1]^2.)/(float(n_ratios-1.)*total(1./abundances[*,1]^2.)))
abund_sigma=max([abund_sigma_gauss,abund_sigma_direct])


print, 'Calculate abundance ratios'
print, 'compared to Allen, Astrophysical quantities (1973)'
print, '-----------------------------------------------------------------------------------------------------------------'
print, 'ratios and error per line'
print, abundances
print, 'MEAN RESULT'
print, format='("abundance ratio: " , e8.2," with error: ", e8.2)',abund,abund_sigma
print, ''
end