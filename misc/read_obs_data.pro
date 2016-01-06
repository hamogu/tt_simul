;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	READ_OBS_DATA
;
; PURPOSE:
;	This procedure reads observational data from a text file
;
; CATEGORY:
;	fitting
;
; CALLING SEQUENCE:
;	READ_OBS_DATA, filename,obsdata,ref
;
; INPUTS:
;	filename:	path and name of the file to be opened, if '' then a dialog opens
;			format: lines with 4 entries (float: photon flux, error, wavelength [angstroem], string:name)
;
; OPTIONAL INPUTS:
;	nh		hydrogen column density in cm^-2 which is used for correction of the observations
; KEYWORDS:
;	no_processing	just reads in the lines of the files in a structure without any modification
;	photons		fluxes are in photons DEFAULT: erg
;
; OUTPUTS:
;	obsdata:	a structure with the following fields:
;			flux	energy flux
;			error	energy error
;			lambda	wavelength in Angstoem, set to zero if energy flux is given
;			name 	a following string used to identify line in data file
;			iz	if name starts with a spectroscopic code -> iz
;			ion	as iz, the ion number
;	ref:		any comments enclosed in '-1' at the end of the file
;
; COMMON BLOCKS:
;	tt_units
;	tt_constants
;
; TEST STATUS:
;	tested
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther 22.09.2005
; 	20.09.07 corrected output bug - lambda was output in cm - MG
;	16.11.07 added iz, ion tags to structure in case the name is interpretable in the first 9 characters
;	03.02.09 added photons keyword
;-


pro read_obs_data, filename,obsdata,ref,nh=nh,no_processing=no_processing,photons=photons
;
common tt_constants
common tt_units

if n_params(0) lt 3 then begin
   print,''
   print,' >   read:obs:data ,filename,data,ref'
   print,''
   return
endif
;
if strtrim(filename,2) eq '' then begin
;
    dir=concat_dir(!xuvtop,'ioneq')
    filename=dialog_pickfile(path=dir,filter='*.dat',title='Select Data File')
    print,' selected:  ',filename
endif
;
;
;
;
openr,lu,filename,/get_lun
;
string1=' '
str=''
;
;
;


while strpos(string1,'-1') EQ  -1 or strpos(string1,'-1') GT 2  do begin
readf,lu,string1
if(strpos(string1,'-1')   EQ  -1 or strpos(string1,'-1') GT 2) then begin
  observed={flux:1.0,error:1.0,lambda:1.0,name:'name'}
  reads,string1,observed
  observed={flux:observed.flux,error:observed.error,lambda:observed.lambda,name:observed.name,iz:0,ion:0}
  obs= n_elements(obs) gt 0 ? [obs,observed]:observed
endif
endwhile
;cut out blanks
;obs.name=strcompress(obs.name,/remove_all)
;apply inverse weight factors, could be e.g. for a transformation photon to energy flux 

if ~keyword_set(no_processing) then begin
  ;transform from photon to energy flux
  obs.lambda*=angstroem ;Angtroem in cm
  index=where(obs.lambda ne 0)
  if ~keyword_set(photons) then obs.flux[index]=obs.flux[index]*c_light*h_planck/obs.lambda[index]
  if ~keyword_set(photons) then obs.error[index]=obs.error[index]*c_light*h_planck/obs.lambda[index]
  for i=0,n_elements(index)-1 do begin
    spectroscopic2ion, strtrim(obs[index[i]].name,1), ionname
    convertname,ionname,iz,ion
    obs[index[i]].iz=iz
    obs[index[i]].ion=ion
  endfor
endif

;
;  get references
refstring=strarr(100)
nref=0
;
string1=' '
while strpos(string1,'-1') EQ  -1  do begin
readf,lu,string1


if(strpos(string1,'-1') EQ -1) and (strpos(string1,'%file') EQ -1) then begin
  refstring(nref)=string1
  nref=nref+1
endif
endwhile
;
ref=refstring(0:nref-1)
;
;
free_lun,lu
;
;
if ~keyword_set(no_processing) then begin
  ;correct for interstellar absorbtion
  obs.lambda/=angstroem ;cm in Angstroem, because output should be in Angstroem
  if n_elements(nh) ne 0 then begin
    tau=exp(-ismtau((findgen(2000)/50),nh=nh,/bam))
    ind=sort(obs.lambda)
    flux=reform((obs.flux)[ind])/spline(findgen(2000)/50,tau,(obs.lambda)[ind],/double)
    obs.flux=flux[sort(ind)]
    error=reform((obs.error)[ind])/spline(findgen(2000)/50,tau,(obs.lambda)[ind],/double)
    obs.error=error[sort(ind)]
  endif

  
endif
obsdata=obs

end
