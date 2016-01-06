;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	line_intens()
;
; PURPOSE:
;	This function calculates the line emissivity 
;	given the different models in the same working directory
;
; CATEGORY:
;	analysis
;
; CALLING SEQUENCE:
;	Result=line_intens, v0,ions,wavelength,abundfile=abundfile,d_lambda=d_lambda,radtemp=radtemp,rphot=rphot,verbose=verbose
;
; INPUTS:
;	v0:	infall velocitie(s), used to identify model files
;	ions:   array of ions of interest in spectroscopic notation
;	wavelength: wavelength in Angstrom for each line, one per ion
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
;	emissivity per volume emission meassure,
;	multiply with n_ion*n_e*V
;	array: ntemp*ndens*nlambda
;
;
; COMMON BLOCKS:
;	elements:	from CHIANTI
;
;
; EXAMPLE:
;		test=line_intens(525.,'C III',977.)
;
; TESTS:
;	TBD
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther, 02.06.2006
;-

function line_intens, v0,ions,wavelength,abundfile=abundfile,d_lambda=d_lambda,radtemp=radtemp,rphot=rphot,verbose=verbose
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

; -- check input for errors ---
if n_elements(v0) eq 0 then begin 
  print,'Specify infall velocities'
  return, -1
endif  
if n_elements(ions) eq 0 then begin
  print,'Specify a string array of ions!'
  return, -1
endif  
if n_elements(ions) ne n_elements(wavelength) then begin
  print,'Number of wavelegth does not match number of ions!'
  return, -1
endif  
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
cd,'./',current=old_dir
descr=[descr,'using models in: '+old_dir]
descr= n_elements(radtemp) eq 0 ? [descr,'No radiative excitation'] : [descr,'Rad Temp: '+string(radtemp)]
descr= n_elements(rphot) eq 0 ? [descr,'No radiative excitation'] : [descr,'dilution factor: '+string(rphot)]
descr=[descr,'Velocity units: km/s']
descr=[descr,'Information about abundance file used']
descr=[descr,abund_ref]
output={infall_vel:v0,lines:lines,description:descr}

; -- place input in structure --
for linenumber=0,n_elements(output.lines)-1 do begin
  output.lines[linenumber].ion=ions[linenumber]
  output.lines[linenumber].wavelength=wavelength[linenumber]
endfor  

for i=0,n_elements(v0)-1 do begin
  print, 'Calculating lines for infall velocity '+string(v0[i])+' km/s'
  restore,strtrim(string(ulong(v0[i])),2)+'.sav'
  read_ioneq,strtrim(string(ulong(v0[i])),2)+'.ioneq',ioneq_logt,ioneq,ioneq_ref
  prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens
  for linenumber=0,n_elements(output.lines)-1 do begin
    if keyword_set(verbose) then print, output.lines[linenumber].ion
    ;decode ion name from sprectroscopic notation to iz, ion as integer numbers
    spectroscopic2ion, output.lines[linenumber].ion, ionname
    convertname,ionname,iz,ion
    intens=calcmygoft(iz,ion,[output.lines[linenumber].wavelength-d_lambda/2.,output.lines[linenumber].wavelength+d_lambda/2.],reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot)
    output.lines[linenumber].intens[i]=total(intens*volem)
  endfor
endfor        

return, output

end

;ions=['He II','C II','C II','C II','C III','C IV','C IV','N III','N III','N V','O VI','O VI','Ne VII','Mg II','Mg II','Si II','Si IV','Fe II','Fe II']
;lamb=[1640.,1334.5,1335.7,2330,977,1548,1550,989.8,991.5,1238,1032,1038,973,2797,2803,1533,1393,2392,2618]
;v0=[250,350,450,525]

;modified add IUE lines, but ignore X II, becaus they form usually at 30.000 K, which is too cold for my model
;ionlist=['C III','C III','C IV','C IV','N III','N III','N V','O VI','O VI','Ne VII','Si III','Si IV','Si IV']
;lambda=[977.,1909,1548,1550,989,991,1238,1032,1038,973,1892,1393,1403]
;


;f=[.25,0,.005,.19,.08,0,.09,.07,.07,0,.000008,.25,.26]
function tt_optdepth, v0,ions,wavelength,f,abundfile=abundfile,d_lambda=d_lambda,radtemp=radtemp,rphot=rphot,verbose=verbose
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

; -- check input for errors ---
if n_elements(v0) eq 0 then begin 
  print,'Specify infall velocities'
  return, -1
endif  
if n_elements(ions) eq 0 then begin
  print,'Specify a string array of ions!'
  return, -1
endif  
if n_elements(ions) ne n_elements(wavelength) then begin
  print,'Number of wavelegth does not match number of ions!'
  return, -1
endif  
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
line_template={ion:'H_I',wavelength:'42.',effective_opt_depth:make_array(n_elements(v0),val=42.)}
lines=replicate(line_template,n_elements(ions))
descr='Line intensities calculated by analysis moduls for TT-Simul'
descr=[descr,'Computed at: '+systime()]
cd,'./',current=old_dir
descr=[descr,'using models in: '+old_dir]
descr= n_elements(radtemp) eq 0 ? [descr,'No radiative excitation'] : [descr,'Rad Temp: '+string(radtemp)]
descr= n_elements(rphot) eq 0 ? [descr,'No radiative excitation'] : [descr,'dilution factor: '+string(rphot)]
descr=[descr,'Velocity units: km/s']
descr=[descr,'Information about abundance file used']
descr=[descr,abund_ref]
output={infall_vel:v0,lines:lines,description:descr}

; -- place input in structure --
for linenumber=0,n_elements(output.lines)-1 do begin
  output.lines[linenumber].ion=ions[linenumber]
  output.lines[linenumber].wavelength=wavelength[linenumber]
endfor  

for i=0,n_elements(v0)-1 do begin
  print, 'Calculating lines for infall velocity '+string(v0[i])+' km/s'
  restore,strtrim(string(ulong(v0[i])),2)+'.sav'
  read_ioneq,strtrim(string(ulong(v0[i])),2)+'.ioneq',ioneq_logt,ioneq,ioneq_ref
  prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens
  for linenumber=0,n_elements(output.lines)-1 do begin
    if keyword_set(verbose) then print, output.lines[linenumber].ion
    ;decode ion name from sprectroscopic notation to iz, ion as integer numbers
    spectroscopic2ion, output.lines[linenumber].ion, ionname
    convertname,ionname,iz,ion
    opt_depth=calculate_opt_depth(wavelength[linenumber],wavelength[linenumber]+[-d_lambda,d_lambda], 200,output.lines[linenumber].ion , hydrodyn,f[linenumber])   
    output.lines[linenumber].effective_opt_depth[i]=opt_depth.effective_opt_depth
  endfor
endfor        

return, output

end
