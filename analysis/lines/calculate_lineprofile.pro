;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	calculate_lineprofile
;
; PURPOSE:
;	as the name says, surprisingly this function is meant to estimate a single lineprofile including doppler effects
;
; CATEGORY:
;	modeling
;
; CALLING SEQUENCE:
;	result=calculate_lineprofile( lambda_0,lambda_range, nbin, snote, hydrodyn,turbulent_broadening=turbulent_broadening,f=f)
;
; INPUTS:
;	lambda_0:	rest wavelength 
;	lambda_range:	2 dim vector containig the lower and upper bound of the wavelength region calculated [Anstroem]
;	If f is not given, lambda_0 and lambda_range can include several lines, lamda_0 is then a vector, lambda_range a 2-dim array
;	nbin:		number of bins
;	snote:		spectroscopic notation of ion (e.g. 'Ne IX')
;	hydrodyn:	from tt_simul
;	
;
; OPTIONAL INPUTS
;	turbulent_broadening:	in Agnstroem, default=0
;	f:		oszillator stregth for specified line
;		If given the optical depth is determeined and a simple I=I_0*exp(-tau) extinction law assumed
;
; OUTPUTS:
;	result: structure with the following tags:
;		.wave		lambda_grid used in Angstroem
;		.flux		flux at stellar surface in erg/cm^2/s/Angstoem
;
; PROGRAMMING NOTES:
;	In priciple the f can be found from the CHIANTI database, give snote and lambda. However this requieres some programming
;	which I felt not necessary as long as only a few lines are analysed this way.
;	The line will be placed as lambda_0, the routine does not check for overlap with other lines.
;
; TEST STATUS:
;	 tested
;
; EXAMPLE:
;	result=calculate_lineprofile(1032.,[1030,1035], 500, 'O VI', hydrodyn,f=2.734e-01)
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 23.11.2005
;	12.02.2007	changed call to calcmygoft to use lambda_range not [lambda_0-1,lambda_0+1], which is no approtiate for all lambda  -HMG
;-
function calculate_lineprofile, lambda_0,lambda_range, nbin, snote, hydrodyn,turbulent_broadening=turbulent_broadening,f=f
;Only one line per interval calculated! Calcmygoft will not resolve different lines
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
;convert from integer
lambda_0=double(lambda_0)
lambda_range=double(lambda_range)

;check input for plausibility
if n_params() lt 5 then begin
   message,"Necessary parameters missing" 
   return,-1
endif

;decode ion name from sprectroscopic notation to iz, ion as integer numbers
spectroscopic2ion, snote, ionname
convertname,ionname,iz,ion
;calculate the intensity for every lambda_0
n_lambda=n_elements(lambda_0)<n_elements(lambda_range)/2
if n_lambda ne 1 then message,'Lineprofile not implemented yet for use with multiples lines',/informational

prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens

;-- core call---
intens=calcmygoft(iz,ion,lambda_range,reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot)
for i=0,n_lambda-1 do intens[*,i]*=volem

;setup output
binwidth=(lambda_range[1,*]-lambda_range[0,*])/nbin
data={wave:make_array(nbin,n_lambda),flux:make_array(nbin,n_lambda)}

;in case extinction shell be considered
if n_elements(f) ne 0 then begin 
  tau=calculate_opt_depth(lambda_0,lambda_range, nbin, snote, hydrodyn,f,turbulent_broadening=turbulent_broadening)
  if n_lambda ne 1 then begin
    message,'Computation for first lambda only.',/continue
    n_lambda=1
  endif
endif    

for i=0,n_lambda-1 do data.wave[*,i]=findgen(nbin)*binwidth[i]+lambda_range[0,i]+0.5*binwidth[i]

for i=0,n_elements(volem)-1 do begin  	;over all layers
  for j=0,n_lambda-1 do begin		;over all lines
    broadening= n_elements(turbulent_broadening) eq 1 ? turbulent_broadening : 0
    broadening +=doppler_broadening(lambda_doppler(lambda_0[j],hydrodyn[3,i]),hydrodyn[4,i],iz)
    if n_elements(f) eq 0 then begin ;optical thin case
      data.flux[*,j]+=intens[i,j]*line_profile_gauss(lambda_doppler(lambda_0[j],hydrodyn[3,i]),data.wave[*,j],broadening);*exp(-1e-13*total(hydrodyn[2,i:*]))
    endif else begin  ;optical depth dims deeper layers
      data.flux[*,j]+=intens[i,j]*line_profile_gauss(lambda_doppler(lambda_0[j],hydrodyn[3,i]),data.wave[*,j],broadening)*exp(-tau.opt_depth[i,*])
    endelse  
  endfor
endfor

return, data
end

