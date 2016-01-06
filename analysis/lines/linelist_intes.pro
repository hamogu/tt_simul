;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	linelist_intens()
;
; PURPOSE:
;	This function calculates the line emissivity 
;	of any number of lines in a list.
;
; CATEGORY:
;	analysis
;
; CALLING SEQUENCE:
;	Result=linelist_intens,ions,wavelength,temp,em,abundfile=abundfile,d_lambda=d_lambda,radtemp=radtemp,rphot=rphot,verbose=verbose
;
; INPUTS:
;	ions:   array of ions of interest in spectroscopic notation
;	wavelength: wavelength in Angstrom for each line, one per ion
;	temp: voctor of temperatures of emission components
;	em: Vector of Emission meassures, default=1
;	
; OPTIONAL INPUTS:
;	abundfile: abundace file to be used default: abund found in common block or !abundfile
;	rphot:	distance from star centre, default is 1=star surface
;	radtemp:temperature of radiation field, default is 6000 K
;	
; KEYWORD PARAMETERS:
;	verbose:	more printout
;
; OUTPUTS:
;	emissivity array: nlambda
;
;
; COMMON BLOCKS:
;	elements:	from CHIANTI
;
;
; EXAMPLE:
;		test=line_intens('C III',977.,3e5)
;
; TESTS:
;	TBD
;
; MODIFICATION HISTORY:
; 	Adapted from line_intes:	Moritz Guenther, 06.09.2007
;-

function line_intens, ions,wavelength,temp, em,abundfile=abundfile,d_lambda=d_lambda,radtemp=radtemp,rphot=rphot,verbose=verbose, dens=dens
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

; -- check input for errors ---
if n_elements(ions) eq 0 then begin
  print,'Specify a string array of ions!'
  return, -1
endif  
if n_elements(ions) ne n_elements(wavelength) then begin
  print,'Number of wavelegth does not match number of ions!'
  return, -1
endif
if n_elements(em) eq 0 then em=temp*0.+1.	;set default value
if n_elements(temp) ne n_elements(EM) then begin
  print,'Array size for temperatures does not match emission meassures!'
  return, -1
endif
if n_elements(em) eq 0 then em=temp*0.+1.e10	;set default value
if keyword_set(abundfile) then begin
  read_abund, abundfile ,abund,abund_ref 
  print,"Using elemental abundances file: "+abundfile
endif  
if n_elements(abund) eq 0 then begin
  read_abund, !abund_file ,abund,abund_ref
  print, "Loading defauld Abundance file: "+!abund_file
endif  

if n_elements(d_lambda) eq 0 then d_lambda=2.


; -- create output structure ---
line_template={ion:'H_I',wavelength:'42.',intens:make_array(n_elements(v0),val=42.),unit:'erg/cm^2/s @ stellar surface'}
lines=replicate(line_template,n_elements(ions))
descr='Line intensities calculated by analysis moduls for TT-Simul'
descr=[descr,'Computed at: '+systime()]
descr= n_elements(radtemp) eq 0 ? [descr,'No radiative excitation'] : [descr,'Rad Temp: '+string(radtemp)]
descr= n_elements(rphot) eq 0 ? [descr,'No radiative excitation'] : [descr,'dilution factor: '+string(rphot)]
descr=[descr,'Information about abundance file used']
descr=[descr,abund_ref]
output={model_em:em, model_temp:temp,lines:lines,description:descr}

; -- place input in structure --
for linenumber=0,n_elements(output.lines)-1 do begin
  output.lines[linenumber].ion=ions[linenumber]
  output.lines[linenumber].wavelength=wavelength[linenumber]
endfor  

for i=0,n_elements(ions)-1 do begin
  for linenumber=0,n_elements(output.lines)-1 do begin
    if keyword_set(verbose) then print, output.lines[linenumber].ion
    ;decode ion name from sprectroscopic notation to iz, ion as integer numbers
    spectroscopic2ion, output.lines[linenumber].ion, ionname
    convertname,ionname,iz,ion
    intens=calcmygoft(iz,ion,[output.lines[linenumber].wavelength-d_lambda/2.,output.lines[linenumber].wavelength+d_lambda/2.],reform(dens),reform(temp),/quiet,radtemp=radtemp,rphot=rphot)
    output.lines[linenumber].intens[i]=total(intens*em)
  endfor
endfor

return, output

end
