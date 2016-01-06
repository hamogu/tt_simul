;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_FB_RAD_LOSS
;
; PURPOSE:
;	This procedure calculates the radiative losses due to free-bound emission and the number of recombinations.
;	
;
;
; CATEGORY:
;	CHIANTI
;
; CALLING SEQUENCE:
; tt_fb_rad_loss, q_fb_rad_loss,t_ion
;
; INPUTS:
;	t_ion:		ion temperature in K
;
; OUTPUTS:
;	q_fb_rad_loss:	radiation loss in erg*cm^3/s/unit emission meassure
;
; TEST STATUS:
;	tested in tt_microscopic and moved to a single routine
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther 29.03.2005
;	Adapted from CHIANTI freebound.pro
;	22.4.05 deleted recomb output
;-
pro tt_fb_rad_loss, q_fb_rad_loss,t_ion,noverner=noverner
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
wvl=double(findgen(1000)+1)
IF n_elements(min_abund) EQ 0 THEN min_abund=0.
nwvl=n_elements(wvl)
int=dblarr(nwvl)
ewvl=1d8/wvl
nt=n_elements(t_ion)

ident_wvl=make_array(nwvl,val=1.)
ident_t=make_array(nt,val=1.)

read_ip,!xuvtop+'/ip/chianti.ip',ionpot,ipref

read_klgfb,pe,klgfb
ksize=size(klgfb)
max_ngf=ksize(2)

IF NOT keyword_set(noverner) THEN BEGIN
  vdata=dblarr(10,465)
  dir=concat_dir(!xuvtop,'continuum')
  fname=concat_dir(dir,'verner_short.txt')
  openr,lun,fname,/get_lun
  readf,lun,vdata
  free_lun,lun
ENDIF
FOR iz=1,n_elements(ioneq[0,*,0]) DO BEGIN		;n_elements(ioneq[*,0]) = number of chemical elements
    ab=abund[iz-1]/total(abund)				;normalize abund
    IF ab Le min_abund THEN GOTO,lbl1
    FOR ion=1,iz DO BEGIN
      ieq=ioneq[0,iz-1,ion]	;get_ieq(temp,iz,ion+1,ioneq_logt=ioneq_logt,ioneq_frac=ioneq)
      ip=ionpot[iz-1,ion-1]
      IF total(ieq) NE 0. THEN BEGIN
        freebound_ion,t_ion,wvl,rad,iz,ion,ip=ip, $
             vdata=vdata,pe=pe,klgfb=klgfb,noverner=noverner
        rad=rad*ab*(ident_wvl#ieq)*1.e-40*4*!pi		;4*pi = emission in all directions
        ;recomb[iz-1,ion]=total(rad*((wvl/1.9865d-8)#ident_t))  ;number of emitted photons=number of recombinations
	int=int+rad
      ENDIF
    ENDFOR
    lbl1:
ENDFOR
q_fb_rad_loss=total(int);*A*delta_x*delta_t*n_e*n_ion	;total(int) is integration dlambda with delta_lamba=1 A
return
end
