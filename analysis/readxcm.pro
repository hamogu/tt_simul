;+
; NAME:
;	readxcm
;
; PURPOSE:
;	This procedure reads XSPEC files which contain fits with phabs(tt_simul_model + vapec +vapec)
;	and write abund files, so the XSPEC results can be used as CHIANTI inputs.
;
; CATEGORY:
;	data conversion
;
; CALLING SEQUENCE:
;	readxcm, xcmfile,abundfile,vapec=vapec,nh=nh
;
; INPUTS:
;	xcmfile:	XSPEC fit file to be read
;
; OUTPUTS:
;	abundfile:	CHIANTI abund file to be written
;
; OPTIONAL OUTPUTS:
;	vapec:		2*2 array with the vapec temp (in K) and normalisation (same untis as XSPEC)
;	nh:		N_H in 1/cm^2
;
; COMMON BLOCKS:
;	tt_elements:	abund and ioneq - from CHIANTI
;
; SIDE EFFECTS:
;	sets !abund_file (CHIANTI) the the basic set of abundances used in XSPEC
;
; TEST STATUS:
;	tested
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 18.08.2007
;-

pro readxcm, xcmfile,abundfile,vapec=vapec,nh=nh
common elements, abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

openr,xcm,xcmfile, /get_lun

dummyline='string'

; --- Parse xcm file ---
; read abund an do some tests
readf,xcm,dummyline	; statistic chi
readf,xcm,dummyline	; abund grsa
dummy=strsplit(dummyline,' ',/extract)
base_abund=dummy[n_elements(dummy)-1]
case base_abund of
  'grsa': defsysv,'!abund_file',!xuvtop + 'abundance/grevesse_sauval98.abund'
  else: message,'Abundance table '+base_abund+' not implemented in readxcm.pro!!!'
endcase
read_abund,!abund_file ,abund,abund_ref
readf,xcm,dummyline	; xsect bcmc
readf,xcm,dummyline	; xset forcecalc off
readf,xcm,dummyline	; cosmo   70.000   0.000   0.730
readf,xcm,dummyline	; model phabs( atable{large_grid.fits} + vapec + vapec )
if ~strmatch(dummyline,'model phabs( atable{*.fits} + vapec + vapec )') then message,'Check that model is phabs* (tt_simul +...)! Otherwise readxcm does not work properly!',/informational

readf,xcm,dummyline	; n_h
reads,dummyline,nh,f1,f2,f3,f4,f5
nh*=1e22		;convert from XSPEC norm to 1/cm^2

readf,xcm,dummyline	; v0
readf,xcm,dummyline	; n0

element=[6,7,8,10,12,14,16,26]  ; which element is in which line of the xcm file
corr=element*0.
for i=0,n_elements(element)-1 do begin
  readf,xcm,dummyline
  if strmid(dummyline,1,1) ne '=' then begin
    reads,dummyline,abund_rel,f1,f2,f3,f4,f5
    abund[element[i]-1]*=abund_rel	
    corr[i]=abund_rel
  endif else begin
    abund[element[i]-1]*= corr[strmid(dummyline,2)-4]
  endelse
endfor
;read temperature and normalisation for the vapec in first data group

readf,xcm,shocknorm,f1,f2,f3,f4,f5
readf,xcm,vapec1kt,f1,f2,f3,f4,f5
for i=0,13 do readf,xcm,dummyline
readf,xcm,vapec1norm,f1,f2,f3,f4,f5
readf,xcm,vapec2kt,f1,f2,f3,f4,f5
for i=0,13 do readf,xcm,dummyline
readf,xcm,vapec2norm,f1,f2,f3,f4,f5
kinkev=8.617e-8
vapec=[[vapec1kt/kinkev,vapec1norm],[vapec2kt/kinkev,vapec2norm]]

free_lun,xcm

; --- Write new abund file ---
openw,af,abundfile,/get_lun
for i=0,29 do if abund[i] gt 0. then printf,af,i+1,alog10(abund[i])+12.
printf,af,'-1'
printf,af,'%file:' ,abundfile
printf,af,'%parsed xcm: ',xcmfile
printf,af,'%comment: produced by readxcm, Moritz Guenther'
printf,af,'-1'
free_lun,af

end