;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	calcmygoft()
;
; PURPOSE:
;	This function calculates the line emissivity per unit mission meassure
;	This supplements mygoft as it loops through all steps. Mygoft works with ioneq array but may fails in locating matching T/ioneq pair if
;	ioneq is not sorted.
;
; CATEGORY:
;	analysis
;
; CALLING SEQUENCE:
;	Result=MYGOFT(IZ,ION,WLOW,WHIGH,LOGDENS,TLOG,quiet=quiet,rphot=rphot,radtemp=radtemp)
;
; INPUTS:
;	iz:	atomic number of elements
;	ion:	spectroscopic number of ion
;	w:	upper and lower wavelength bound in Angstrom (w[0,*] is lower bound and w[1,*] upper)
;	dens:	n_e [1/cm^3]
;	Temp:	Temperature in Kelvin
;
; OPTIONAL INPUTS:
;	rphot:	distance from star centre, default is 1=star surface
;	radtemp:temperature of radiation field, default is 6000 K
;	
; KEYWORD PARAMETERS:
;	quite:	suppresses printout
;
; OUTPUTS:
;	emissivity per volume emission meassure,
;	multiply with n_ion*n_e*V
;
;
; COMMON BLOCKS:
;	elements:	from CHIANTI
;
; SIDE EFFECTS:
;	sets rphot and radtemp to standard values if undefined previously
;
; TESTS:
;	
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 17-05-2005
;			took factors of ioneq, etc out of routine, implemented rphot and radtemp
;-
function calcmygoft, IZ,ION,W,DENS,Temp,quiet=quiet,rphot=rphot,radtemp=radtemp
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
;backup ioneq etc
fullioneq=ioneq
if n_elements(ioneq_logt) gt 0 then full_ioneq_logt=ioneq_logt
nd=n_elements(dens)
nt=n_elements(temp)
nlambda=n_elements(w[0,*])
if (nd eq nt and n_elements(ioneq[*,0,0]) eq nd) then begin
  ;fill result with negative values, so errors are easier to detect
  result=make_array(nd,nlambda,val=-1.)
  for i=0,n_elements(temp)-1 do begin
    ;set ioneq to a sinlge value to prevend unwanted interpolations
    ioneq_logt=alog10(temp[i])
    ioneq=fullioneq[i,*,*]
    result[i,*]=mygoft(IZ,ION,W,alog10(dens[i]),alog10(temp[i]),quiet=quiet,rphot=rphot,radtemp=radtemp)
  endfor
  ioneq=fullioneq
  if n_elements(full_ioneq_logt) gt 0 then ioneq_logt=full_ioneq_logt
  return, result
endif else begin
  message,'ioneq, Density and temperature are not of the same size'
  return, -1
endelse  
end