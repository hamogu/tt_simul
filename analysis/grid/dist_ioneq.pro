;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	dist_ioneq()
;
; PURPOSE:
;	This function calculates a number which is meassure for the deviation of the ionisation state from equilibrium
;
; CATEGORY:
;	analysis
;
; CALLING SEQUENCE:
;	Result=dist_ioneq(ionstate,t_e,[ioneq_file=ioneq_file, static=static])
;
; INPUTS:
;	ionstate:	the ionstate in question in CHIANTI ioneq format
;	T_e:		electron temperature in K
;
; OPTIONAL INPUTS:
;	ioneq_file:	a file with equilibrium ionstates, default: calculate equilib with tt_test_ioneq
;	static:		equilibrium ioneq, saves time if dist_ioneq is called multiple times and static ia already defined
;
; KEYWORDS:
;	elements:	output is an 2-dim array, the second dim is distace per element
;	
; OUTPUTS:
;	a number between 0 and 1. 0 means a complete match with the equilibrium, 1 that the is no ionisation state that is populated
;	in the this case and the equilibrium. It is calculated as the absolute deviation from equilibrium and averaged ov3er all elements.
;	The abundance is not taken into account. One could argue that elements should be weighted by their abundance but this would lead to 
;	a result completly dominated by hydrogen. For some questions that is OK but for others it is not. For the lines in the X-ray region
;	C,O, Ne or Fe are more impotant. In oder the represent the deviation from equilibrium by one number only this definition was chosen.
;
; EXAMPLE:
;		IDL> setup_elements
;		IDL> dist_ioneq(ioneq, 6e3)
;	Should deliver a small number because data loaded by setup_elements shell represent equilibria
;
; TESTS:
;	tested visually
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther, 11.7.05
;	5.10.05		added /elements
;-
function dist_ioneq, ionstate,t_e,ioneq_file=ioneq_file, static=static,elements=elements
if n_elements(ionstate[*,0,0]) ne n_elements(t_e) then begin
  print, 'T_e must contain one temperature value for each ion state!'
  return, -1
endif else begin 
  
  if n_elements(ioneq_file) eq 1 then begin
    read_ioneq,ioneq_file ,ioneq_logt,ioneq_frac,ioneq_ref
  endif else begin
    if n_elements(static) eq 0 then tt_test_ioneq,static,/quiet
    ioneq_logt=findgen(41)/10.+4.
    ;due to small numerical errors in lsode there may be negative values in static of the order -1e-8
    static=abs(static)
  endelse
  n_elem=n_elements(ionstate[0,*,0]) 
  dist=make_array(n_elements(t_e),n_elem,val=0.)
  for z=0,n_elem-1 do begin
    for ion=0,z+1 do begin
      ;Should this be weighted by abundance?
      dist[*,z]+=abs(get_ieq(t_e, z+1,ion+1,ioneq_logt=ioneq_logt,ioneq_frac=static)-ionstate[*,z,ion])
      ;print, dist[140], ion,z
    endfor
  endfor
  if ~keyword_set(elements) then dist=total(dist,2)
  return, dist
endelse
end      