;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	ANALYSGRID_extra
;
; PURPOSE:
;	This procedure calculates the itensity of a single line from model files 
;
; CATEGORY:
;	analysis
;
; CALLING SEQUENCE:
;	analysegrid_extra,namelist,data,iz,ion,wvl,radtemp=radtemp,rphot=rphot,no_fe=no_fe,filename=filename,no_n=no_n
;
; INPUTS:
;	namelist:	This procedure reads model files as input. namelist is an array with the path to all .sav files which shell be used 
;	iz:		Element number -1 (H is 0, He is 1,...)
;	ion:		charge of ion in question
;
; OPTIONAL INPUTS:
;	radtemp:	temperatue of a background black-body radiation field [Kelvin]
;	rphot:		dilution factor, describing the dilutes of the background radiation field [1 is stellar surface]
;	filename:	if set the data structur is saved to filename after every step
;
; OUTPUTS:
;	data:		a strucure with different fields which contains the results, exspecially the line emissivity in erg/s/cm^2
;			use help, data, /structure for a list of fields
;
; TEST STATUS:
;	checked 
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 22.09.2005
;-
pro analysegrid_extra,namelist,data,iz,ion,wvl,radtemp=radtemp,rphot=rphot
num=n_elements(namelist)
data={v0:fltarr(num),n0:fltarr(num),line_extra: fltarr(num)  }
;namelist is an array of positions for *.sav files
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
read_abund, !abund_file ,abund,abund_ref
    
for i=0,num-1 do begin
  print, 'Working on '+namelist[i]
  restore,namelist[i]
  read_ioneq, strmid(namelist[i],0,strlen(namelist[i])-4)+'.ioneq' ,ioneq_logt,ioneq,ioneq_ref
  prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens
  n=n_elements(hydrodyn[0,*])-1
   ;--- starting hydrodynamics ---
  data.v0[i]		=hydrodyn[3,n]
  data.n0[i]		=hydrodyn[2,n]
  
  data.line_extra[i]=total(volem*calcmygoft(iz,ion,wvl,reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot))
  
  if n_elements(filename) ne 0 then begin 
    save, data,filename=filename
    print, systime()
  endif  
endfor          
end             