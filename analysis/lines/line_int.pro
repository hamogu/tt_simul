;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	line_int
;
; PURPOSE:
;	shows line intensities given a grid of models with different filling facor, e.g. on an inhomogenous spot
;
; CATEGORY:
;	spectral synthsis
;
; CALLING SEQUENCE:
;	line_int,data,snote,wvl,f,distance,quiet=quiet,intensities=intensities,d_lambda=d_lambda
;
; INPUTS:
;	data	:structure as output by spectral_series
;	snote:	spectroscopic notation of the ions in question
;	wvl:	wavelength in Angstroem
;
; OPTIONAL INPUTS:
;	f:	array of filling factors
;	distance: distance to star in pc
; KEYWORDS:
;	d_lambda:  delta_lambda (true wavelength -- input wavelength), default=1.
;
; OPTIONAL OUTPUTS:
;	intensities: 	array of line intensities (n_vel*n_lines)
;	model:		array of line intensities folded with filling factors (n_lines)
;
; COMMON BLOCKS:
;	elements:	abundance and ioneq
;
; EXAMPLE:
;	line_int,uv,['C III','O VI','O VI','N V','C IV','C IV','Si IV','Si IV','S VI','C III'],[977,1032,1038,1238,1548,1550,1394,1403,933,1176],[.002,.01,.05,.15]
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

;line_int,uv,['C III','O VI','O VI','N V','C IV','C IV'],[977,1032,1038,1238,1548,1550],[.002,.01,.05,.15]
pro line_int,data,snote,wvl,f,distance,quiet=quiet,intensities=intensities,d_lambda=d_lambda,model=model
common tt_units
if n_elements(snote) ne n_elements(wvl) then begin
  print, '%spectral_series::line_int : snote and wvl do not have same number of elements. Aborting.'
  return
endif
n_vel=n_elements(data.infall_vel)
n_lines=n_elements(snote)
intensities=make_array(n_vel,n_lines)

case n_params() of
  3: begin
     scaling=1.
     f=make_array(n_vel,val=1.)
     end
  4: begin
     scaling=1.
     end
  5: begin
     ;from distance in pc to scaling factor
     distance=double(distance)
     scaling=(7e10*cm/(distance*pc))^2.
     end
  else: begin
     print,'%spectral_series::line_int : Wrong number of arguments'
     return
     end
endcase
     

for i=0,n_vel-1 do begin
  for j=0,n_lines-1 do begin
    index=get_line_index((*data.scenario[i].ch_spec).lines,snote[j],wvl[j],d_lambda=d_lambda)
    intensities[i,j]=index[0] lt 0? -1: total((*data.scenario[i].ch_spec).lines[index].int)		;index[0] not just index because cannot compare vector to 0
  endfor
endfor
;intensities is now per unit angle. To get total emission multiply by 4 pi
intensities*=4*!pi*scaling
;if f is set, calculate total models
if n_params() ge 4 then begin
  model=make_array(n_lines)
  for j=0,n_lines-1 do model[j]=total(intensities[*,j]*f)
endif  
  
;intensities can be read out as a keyword, but mostly this procedure will be used to print results on screen
if ~keyword_set(quiet) then begin
  print,'Line intensities' 
  if scaling eq 1 then print,'Intensities @ stellar surface' else print,distance,format='("Intensities at distance:",(g6.3)," pc")'
  print,snote,format='("v_0   n_0      |   f   |",'+string(n_lines)+'((A7)," "))'
  print,wvl,format='(  "               |       |"'+string(n_lines)+'((F7.1)," "))'
  print,               "km/s  cm^{-3}  |       |     -------------- erg/(s cm^2) -------------"
  for i=0,n_vel-1 do print, data.infall_vel[i],min((*data.scenario[i].ch_spec).model_ne),f[i],intensities[i,*],format='((f5.0)," ",(e8.1)," | ",(d5.3)," | ",'+string(n_lines)+'((e7.1)," "))'
  if n_params() ge 4 then print, model,format='("f - model      | sum   | "'+string(n_lines)+'((e7.1)," "))'
endif

end  

