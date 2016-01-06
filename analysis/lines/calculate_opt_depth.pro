;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	calculate_opt_depth
;
; PURPOSE:
;	as the name says, surprisingly this function is meant to estimate the optical depth for any given line
;
; CATEGORY:
;	modeling
;
; CALLING SEQUENCE:
;	result=calculate_opt_depth( lambda_0,lambda_range, nbin, snote, hydrodyn,f,turbulent_broadening=turbulent_broadening)
;
; INPUTS:
;	lambda_0:	rest wavelength 
;	lambda_range:	2 dim vector containig the lower and upper bound of the wavelength region calculated [Anstroem]
;	nbin:		number of bins
;	snote:		spectroscopic notation of ion (e.g. 'Ne IX')
;	hydrodyn:	from tt_simul
;	f:		oszillator stregth for specified line
;
; OPTIONAL INPUTS
;	turbulent_broadening:	in Agnstroem, default=0
;
; OUTPUTS:
;	result: structure with the following tags:
;		The first and second are meant for plotting e.g. IDL>plot, result.depth,result.opt_simple_depth
;		.depth		shock front=0, meassured inwards in [cm]
;		.wave		lambda_grid used in Angstroem
;
;		.opt_simple_depth: optical depth tau(z) z, corresponds to .depth
;				This just adds up the optical depth at line centre.
;		.opt_depth:	Because of doppler shift the line centre changes ots lambda. Considering this means 
;				the wings of the gauss profile have to be modeled, too. output: array[depth*wavelegth  grid]
;
; PROGRAMMING NOTES:
;	In priciple the f can be found from the CHIANTI database, given snote and lambda. However this requires some programming
;	which I felt not necessary as long as only a few lines are analysed in this way.
;	The line will be placed as lambda_0, the routine does not check for overlap with another line.
;
; TEST STATUS:
;	 tested
;
; EXAMPLE:
;	result=calculate_opt_depth( 1032.,[1030,1035], 500, 'O VI', hydrodyn,.1367)
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 23.11.2005
;	14.06.2006 Moritz Guenther added mean opt depth to output
;	13.02.2006 Moritz Guenther corrected sqrt(pi) to pi in opt_depth[i,*]=
;-
function calculate_opt_depth, lambda_0,lambda_range, nbin, snote, hydrodyn,f,turbulent_broadening=turbulent_broadening
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
common tt_constants
;--Input plausibility checks
If n_elements(abund) eq 0 or n_elements(ioneq) eq 0 then begin
   message,"Common block elements not initialised" 
   return,-1
endif
If n_elements(ioneq[*,0,0]) ne n_elements(hydrodyn[0,*]) then begin
   message,"Ioneq and hydrodyn not compatible" 
   return,-1
endif
;convert from integer
lambda_0=double(lambda_0)
lambda_range=double(lambda_range)

;check input for plausibility
if n_params() lt 6 then begin
   message,"Necessary parameters missing" 
   return,-1
endif

n_steps=n_elements(hydrodyn[0,*])

;from spectroscopic notation to ion,iz
spectroscopic2ion, snote, ionname
convertname,ionname,iz,ion

prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens
;setup results and wavegrid
binwidth=(lambda_range[1]-lambda_range[0])/nbin
wavegrid=findgen(nbin)*binwidth+lambda_range[0]+0.5*binwidth
result={opt_depth:make_array(n_steps,nbin),opt_simple_depth:make_array(n_steps),depth:make_array(n_steps),wave:make_array(nbin),effective_opt_depth:0.}
result.depth=hydrodyn[0,*]
binwidth=(lambda_range[1,*]-lambda_range[0,*])/nbin
result.wave=findgen(nbin)*binwidth[0]+lambda_range[0]+0.5*binwidth[0]
;temporary variable for the extinction per bin
simple_depth=make_array(n_steps)
opt_depth=make_array(n_steps,nbin)

;start at hock and work inwards - this allows to add up the optical depth for deeper layers immediatly
for i=n_steps-1,0,-1 do begin
  ;optical depth at line centre
  broadening= n_elements(turbulent_broadening) eq 1 ? turbulent_broadening : 0
  broadening +=doppler_broadening(lambda_doppler(lambda_0,hydrodyn[3,i]),hydrodyn[4,i],iz)
  test=min(abs(wavegrid-lambda_0),min_index)
  simple_depth[i]=2.82e-13*!pi*(lambda_0^2*1e-8)*(line_profile_gauss(lambda_0,wavegrid,broadening))[min_index]*f*abund[iz-1]*ioneq[i,iz-1,ion-1]*hydrodyn[2,i]*delta_x[i]
  ;more detailed: lambda_0=lambda_0(z)
  lambda_0_doppler=lambda_doppler(lambda_0,hydrodyn[3,i])
  broadening= n_elements(turbulent_broadening) eq 1 ? turbulent_broadening : 0
  broadening +=doppler_broadening(lambda_doppler(lambda_0_doppler,hydrodyn[3,i]),hydrodyn[4,i],iz)
  opt_depth[i,*]=2.82e-13*!pi*(lambda_0_doppler^2*1e-8)*(line_profile_gauss(lambda_0_doppler,wavegrid,broadening))*f*abund[iz-1]*ioneq[i,iz-1,ion-1]*hydrodyn[2,i]*delta_x[i]
  ; addup the extiction from individual bins to get the optical depth
  result.opt_simple_depth[i]=total(simple_depth[i:*])
  result.opt_depth[i,*]=total(opt_depth[i:*,*],1)
  
  result.effective_opt_depth+=total(result.opt_depth[i,*]*(line_profile_gauss(lambda_0_doppler,wavegrid,broadening))/total(line_profile_gauss(lambda_0_doppler,wavegrid,broadening))*volem[i]*ioneq[i,iz-1,ion-1])/(total(volem*ioneq[*,iz-1,ion-1]))
  
endfor  

; addup the extiction from individual bins to get the optical depth 
;for i=1, n_steps do begin
;  result.opt_simple_depth[n_steps-i]=total(simple_depth[n_steps-i:*])
;  result.opt_depth[n_steps-i,*]=total(opt_depth[n_steps-i:*,*],1)
;endfor  

;get EM weigthed mean opt depth

return, result
end
