;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	mygoft()
;
; PURPOSE:
;	This function calculates the line emissivity per unit mission meassure
;
; CATEGORY:
;	analysis
;
; CALLING SEQUENCE:
;	Result=MYGOFT,IZ,ION,WLOW,WHIGH,LOGDENS,TLOG,quiet=quiet,rphot=rphot,radtemp=radtemp
;
; INPUTS:
;	iz:	atomic number of elements
;	ion:	spectroscopic number of ion
;	w: 	bound of wavelength in Angstrom (first row is lower, second row higher bound)
;	logdens:log10(n_e)
;	tlog	log10(Temperature in Kelvin)
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
;	array: ntemp*ndens*nlambda
;
;
; COMMON BLOCKS:
;	elements:	from CHIANTI
;
; SIDE EFFECTS:
;	sets rphot and radtemp to standat values if undifined previously
;
;
; EXAMPLE:
;		result=mygoft(8,7,[22.,22.2],11.,6.3)
;
; TESTS:
;	can reproduce Fig. 9 from Ness et al. A&A 387,1032-1046(2002)
;	testo7f=mygoft(8,7,22.,22.2,findgen(30)/10.+9.,ioneq_logt,rphot=0.,radtemp=0.)
;	testo7i=mygoft(8,7,21.7,21.9,findgen(30)/10.+9.,ioneq_logt,rphot=0.,radtemp=0.)
;	plot f to i
;
; MODIFICATION HISTORY:
; 	Written by:	Juergen Schmitt
;	13-05-2005	heavyly modified by Moritz Guenther
;			implemented rphot and radtemp
;-
FUNCTION MYGOFT,IZ,ION,W,LOGDENS,TLOG,quiet=quiet,rphot=rphot,radtemp=radtemp

common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

abund_curr=abund(iz-1)

ioneq_ion=REFORM(ioneq(*,iz-1,ion-1))
ind=WHERE(ioneq_ion NE 0.)

if n_elements(rphot) eq 0 then rphot = 1.
if n_elements(radtemp) eq 0 then radtemp = 6000.

ntemp = n_elements(tlog)
ndens = n_elements(logdens)
nlambda=n_elements(w[0,*])
cal_emiss = fltarr(ntemp,ndens,nlambda)
cal_emiss[*] = 0.

;emiss_calc returns a structure with lots of information abourt the line emission
repeat begin
  emiss=emiss_calc(iz,ion,temp=tlog,dens=logdens,noprot=noprot,radtemp=radtemp,rphot=rphot,/no_setup,quiet=quiet) 
  test_emiss=size(emiss)
endrep until test_emiss[n_elements(test_emiss)-2] eq 8
;sometimes there are errors in pop_solver (splstr must be a structure in this context). I cannot find the origin as splstr is part of some common block an I have no idea
;where it is actually set. This seems to happen more when using remote PC's and often during backup time. So I suppose it has something to do with network response time or something. As a workaroud I put a check in emiss_calc to detect that and return here where the calculation is just restarted. So I loose only some time, but the program does not stop completly.

for i=0, nlambda-1 do begin
  wlow=w[0,i]
  whigh=w[1,i]
;take only lines between wlow and whigh
  gd = where(emiss.lambda gt WLOW and emiss.lambda lt WHIGH,ct)
  emissw = emiss[gd]
  em    = emissw.em
;if logdens or tlog are not arrays but scalar em will have the wrong number of dimensions
  em=reform(em,ntemp,ndens,n_elements(em)/(ntemp*ndens))

  if not keyword_set(quiet) then print,'Number of lines: ',ct

  for j=0,ntemp-1 do begin		;loop over all temperatures in grid
    for k=0,ndens-1 do begin		;loop over all densities in grid
       cal_emiss[j,k,i] = total(em[j,k,*])	;sum over all levels of ion
    endfor
  endfor
endfor  
func = cal_emiss
for  k=0,ndens-1 do func[*,k,*]=cal_emiss[*,k,*]/10.^(logdens[k])*ioneq_ion[0] * abund_curr ;func(*,k)=0.83*cal_emiss(*,k)*ioneq/10.^(logdens(k)) * abund_curr
zero = where(func le 0.,ct)
if (ct gt 0) then func[zero] = 0.
return,reform(func)
;we return the line emissivity in units of erg/sec per unit volume emission
;measure for that abundance
end