;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	ANALYSGRID
;
; PURPOSE:
;	This procedure calculates the hydrodynamic and spectroscopic key data from a list of model files. 
;
; CATEGORY:
;	analysis
;
; CALLING SEQUENCE:
;	analysegrid,namelist,data,radtemp=radtemp,rphot=rphot,no_fe=no_fe,filename=filename,no_n=no_n
;
; INPUTS:
;	namelist:	This procedure reads model files as input. namelist is an array with the path to all .sav files which shell be used 
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
; KEYWORDS:
;	no_fe:		calculate no iron lines to save time
;	no_n:		calculate no nitrogen lines to save time
;
; TEST STATUS:
;	checked 
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 22.09.2005
;-
pro analysegrid,namelist,data,radtemp=radtemp,rphot=rphot,no_fe=no_fe,filename=filename,no_n=no_n
description=['v0','n0','Te_0','Tion_0','maxdepth','Te_max','Tion_max','ne9f2i','o7f2i','o7f','o7i','o7r','ne9f','ne9i','ne9r']
num=n_elements(namelist)
data={v0:fltarr(num),n0:fltarr(num),Te_0:fltarr(num),Tion_0:fltarr(num),maxdepth:fltarr(num),Te_max:fltarr(num),Tion_max:fltarr(num),$
    ne9f2i:fltarr(num),o7f2i:fltarr(num),o7f:fltarr(num),o7i:fltarr(num),o7r:fltarr(num),ne9f:fltarr(num),ne9i:fltarr(num),ne9r:fltarr(num),$
    olya:fltarr(num),olyb:fltarr(num),nelya:fltarr(num),nlya:fltarr(num),clya:fltarr(num),n_steps:fltarr(num),diff2eq:fltarr(num),$
    fe17l1501:fltarr(num),$
    fe17l1705:fltarr(num),$
    fe17l1709:fltarr(num),$
    fe17l1678:fltarr(num),$
    fe17l1526:fltarr(num),$
    fe18l1420:fltarr(num),$
    fe18l1600:fltarr(num),$
    fe18l1607:fltarr(num),$
    fe18l1467:fltarr(num),$
    fe19l1352:fltarr(num),$
    fe19l1380:fltarr(num),$
    fe19l1466:fltarr(num),$
    n6r:fltarr(num),$    
    n6i:fltarr(num),$    
    n6f:fltarr(num)    }
;namelist is an array of positions for *.sav files
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
read_abund, !abund_file ,abund,abund_ref
tt_test_ioneq,static,/quiet
    
for i=0,num-1 do begin
  print, 'Working on '+namelist[i]
  restore,namelist[i]
  read_ioneq, strmid(namelist[i],0,strlen(namelist[i])-4)+'.ioneq' ,ioneq_logt,ioneq,ioneq_ref
  prepare_analysis, hydrodyn, delta_x,volem,elecdens,o7i_p,o7f_p,o7r_p,o7f2i_p,ne9i_p,ne9f_p,ne9r_p,ne9f2i_p,/quiet,radtemp=radtemp,rphot=rphot
  n=n_elements(hydrodyn[0,*])-1
  
  ;--- starting hydrodynamics ---
  data.v0[i]		=hydrodyn[3,n]
  data.n0[i]		=hydrodyn[2,n]
  data.te_0[i]		=hydrodyn[5,n]
  data.tion_0[i]	=hydrodyn[4,n]
  ;--- shock properties ---
  data.maxdepth[i]	=hydrodyn[0,0]
  data.te_max[i]	=max(hydrodyn[5,*])
  data.tion_max[i]	=max(hydrodyn[4,*])
  data.n_steps[i]	=n_elements(hydrodyn[0,*])
  data.diff2eq[i]	=total(dist_ioneq(ioneq,reform(hydrodyn[5,*]),static=static)/30./2*delta_x)/data.maxdepth[i]
  ;--- line ratios ---
  data.ne9f2i[i]	=total(ne9f_p*volem)/total(ne9i_p*volem)
  data.o7f2i[i]		=total(o7f_p*volem)/total(o7i_p*volem)
  ;--- line strenghts ---
  data.o7f[i]	=total(o7f_p*volem)
  data.o7i[i]	=total(o7i_p*volem)
  data.o7r[i]	=total(o7r_p*volem)
  data.ne9f[i]	=total(ne9f_p*volem)
  data.ne9i[i]	=total(ne9i_p*volem)
  data.ne9r[i]	=total(ne9r_p*volem)
  oly		=calcmygoft(8,8,[[18.96,18.98],[16.,16.1]],reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot)
  data.olya[i]	=total(oly[*,0]*volem)
  data.olyb[i]	=total(oly[*,1]*volem)
  data.nelya[i]	=total(volem*calcmygoft(10,10,[12.13,12.14],reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot))
  if ~keyword_set(no_n) then begin
    data.nlya[i]=total(volem*calcmygoft(7,7,[24.77,24.8],reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot))
    n6		=calcmygoft(7,6,[[28.7,28.9],[29.0,29.2],[29.5,29.6]],reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot)
    data.n6r[i]	=total(n6[*,0]*volem)
    data.n6i[i]	=total(n6[*,1]*volem)
    data.n6f[i]	=total(n6[*,2]*volem)
  endif
  data.clya[i]	=total(volem*calcmygoft(6,6,[33.7,33.75],reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot))
  
  ;   -iron-
  if ~keyword_set(no_fe) then begin
    fe17list	=[15.015,17.0532,17.0981,16.7778,15.2621]
    fe18list	=[14.205,16.005,16.072,14.6705] ;14.205 is for two lines (14,203 and 14,2078)
    fe19list	=[13.5207,13.8001,14.668]	;13.5207 is for [13.5206,13.5208] and 13.8001 for the close lines[13.7937,13.7950,13.8001,13.8096] and 14.668 for [14.6671,14.668,14.6711]
    fe17	=calcmygoft(26,17,transpose([[fe17list-.02],[fe17list+.02]]),reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot)
    fe18	=calcmygoft(26,18,transpose([[fe18list-.01],[fe18list+.01]]),reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot)
    fe19	=calcmygoft(26,19,transpose([[fe19list-.01],[fe19list+.01]]),reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot)
    data.fe17l1501[i]=total(fe17[*,0]*volem)
    data.fe17l1705[i]=total(fe17[*,1]*volem)
    data.fe17l1709[i]=total(fe17[*,2]*volem)
    data.fe17l1678[i]=total(fe17[*,3]*volem)
    data.fe17l1526[i]=total(fe17[*,4]*volem)
    data.fe18l1420[i]=total(fe18[*,0]*volem)
    data.fe18l1600[i]=total(fe18[*,1]*volem)
    data.fe18l1607[i]=total(fe18[*,2]*volem)
    data.fe18l1467[i]=total(fe18[*,3]*volem)
    data.fe19l1352[i]=total(fe19[*,0]*volem)
    data.fe19l1380[i]=total(fe19[*,1]*volem)
    data.fe19l1466[i]=total(fe19[*,2]*volem)
  endif
  if n_elements(filename) ne 0 then begin 
    save, data,filename=filename
    print, systime()
  endif  
endfor          
end             