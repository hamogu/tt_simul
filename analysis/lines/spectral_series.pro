;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	spectral_series
;
; PURPOSE:
;	Calculates spectra using the ch_synthetic engine for a group of simulations, e.g. different infall velocities of the same spot
;
; CATEGORY:
;	spectral synthsis
;
; CALLING SEQUENCE:
;	Result=spectral_series(v0,wavelength,abundfile=abundfile,d_lambda=d_lambda,radtemp=radtemp,rphot=rphot)
;
; INPUTS:
;	v0	:array of infall velocities, which is also usedto identifiy files, e.g. 250.ioneq and 250.hydrodyn
;	wavelength: two element array. First element is lower wavelegth bound, second is upper bound
;
; OPTIONAL INPUTS:
;	abundfile: abundance file in CHIANTI format to be used, default !abund_file
;	d_lambda:  bin width in Angstroem
;	radtemp:   temperature of radiation background field
;	rphot:     distance to stellar surface giving the dilution for the background radiation field
;
; KEYWORDS:
;	continuum:	add continuum spectrum, time intesive!
;
; OUTPUTS:
;	a complicated structure
;	.infall_vel:	array of invall velocities
;	.wave:		array of wavelength in Agstroem of simulated spectra
;	.description:	text array with further information
;	.scenario:	array of structures
;		.v0:	infall vel of this scenario
;		.CH_SPEC:pointer to a structure as it is output from CH_synthetic
;
; COMMON BLOCKS:
;	elements:	abundance and ioneq
;
; EXAMPLE:
;	testuv=spectral_series([525,450,350,250],[900.,1200.],d_lambda=.1,abund='../../tt_data/twhya2.abund')
;
; NOTES:
;	maybe correct pops needs to be replaced by a proxy, reducing accuracy because it does not work with the TT_simul was of treating an ioneq
;
; TEST STATUS:
;	tested
;
; MODIFICATION HISTORY:
;	Written by:	Moritz Günther 07.04.2006
;-


;testuv=spectral_series([525,450,350,250],[900.,1200.],d_lambda=.1,abund='../../tt_data/twhya2.abund')
;eine Linie Fe XVIII bei 975 ist für 525 und 450 Größenordnungen über dem Ergebnis von line_intens, sie ist auch nicht beobachtet. Ursache unbekannt. Igendwie
;Probleme mit erstem oder letzem bin (aber da ioneq quasi 0)? Könnte andere Linien geben, deren intens os niedrig ist, dass ich das Problem da nciht bemerkt habe.
;Das ist aber reproduzierbar, ich ahbe das noch mal laufen lassen mit
;testfe=spectral_series([525,450,350,250],[974.5,975],d_lambda=.1,abund='../../tt_data/twhya2.abund')
;da reingehen und malk schauen warum der Unterschied passiert
; test=line_intens([250,350,450,525],['C III','Fe XVIII'],[977.,975],abund='../../tt_data/twhya2.abund');
; Lösung : correct_pops arbeitet nicht richtig mit dem Weg wie ich mit dem ioneq umgehe

function spectral_series,v0,wavelength,abundfile=abundfile,d_lambda=d_lambda,radtemp=radtemp,rphot=rphot,continuum=continuum
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

; -- check input for errors ---
if n_elements(v0) eq 0 then begin 
  print,'Specify infall velocities'
  return, -1
endif  
if n_elements(wavelength) ne 2 then begin
  print,'Specify [lower wavelegth, upper wavelength] !'
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

if n_elements(d_lambda) eq 0 then d_lambda=float(wavelength[1]-wavelegth[0])/500


; -- create output structure ---
wavegrid=wavelength[0]+findgen((wavelength[1]-wavelength[0])/d_lambda)*d_lambda

scenario={v0:42.,ch_spec:ptr_new(/allocate_heap)}
scenarios=replicate(scenario,n_elements(v0))
;replicate replicates the pointers, so they all point to the same object
for i=1,n_elements(v0)-1 do scenarios[i].ch_spec=ptr_new(/allocate_heap)
descr='Line intensities calculated by analysis moduls for TT-Simul'
descr=[descr,'Computed at: '+systime()]
cd,'./',current=old_dir
descr=[descr,'using models in: '+old_dir]
descr= n_elements(radtemp) eq 0 ? [descr,'No radiative excitation'] : [descr,'Rad Temp: '+string(radtemp)]
descr= n_elements(rphot) eq 0 ? [descr,'No radiative excitation'] : [descr,'dilution factor: '+string(rphot)]
descr=[descr,'Velocity units: km/s']
descr=[descr,'Information about abundance file used']
descr=[descr,abund_ref]
output={infall_vel:v0,wave:wavegrid,scenario:scenarios,description:descr}

; -- place input in structure --
for scenarionumber=0,n_elements(v0)-1 do begin
  output.scenario[scenarionumber].v0=v0[scenarionumber]
endfor  

; -- place input in structure --
for i=0,n_elements(v0)-1 do begin
  print, 'Calculating lines for infall velocity '+string(v0[i])+' km/s'
  restore,strtrim(string(ulong(v0[i])),2)+'.sav'
  read_ioneq,strtrim(string(ulong(v0[i])),2)+'.ioneq',ioneq_logt,ioneq,ioneq_ref
  prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens
  tt_synthetic, wavelength[0],wavelength[1],out=temp,dens=hydrodyn[2,*],logt_i=ioneq_logt,ioneq=strtrim(string(ulong(v0[i])),2)+'.ioneq',/tt_sim,logem_i=alog10(volem)
  ;no need to store tt_syn, as all information is kept in tempspec anyway
  make_chianti_spec, temp,wavegrid,tempspec,instr=0.06,abund=abundfile,continuum=continuum
  *output.scenario[i].ch_spec=tempspec
endfor        

return, output

end