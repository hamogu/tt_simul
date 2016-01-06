;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_TEST_IONEQ_RATE
;
; PURPOSE:
;	test tt_microscopic
;	For constant temperature and long times tt_microscopic should convert towards 
;	an ionisation equilibrium - in this case Mazzotta et al, where my rates come from
;
; CATEGORY:
;	VERIFICATIOM
;
; CALLING SEQUENCE:
;	TT_test_ioneq_rate,ion_file=ion_file,outpath, abundfile=abundfile, temperature=temperature,result=result
;      
; INPUTS:
;	outpath:	path to write data and log files with prefix, eg. './test25'
;
; OPTIONAL INPUTS:
;	abundfile:	abundance file, if not specified !abund_file is used
;	temperature:	If specified, test is only done for a single temperature
;	ion_file:	The .ioneq file to be used as starting point (in CHIANTI format)
;			The result should not depend on this as long as all elements are included.
;			If not set test5.ioneq is used which contains all elements
;      
; OPTIONAL OUTPUTS:
;	result:		an *.ioneq array as result of simulation (without /temperature)
;			or a time series to check convergence (with temperature=...)
;                        
; EXAMPLE:
;	tt_test_ioneq_rate,'../tt_results/test24'
;
; TEST STATUS:  
;	tested and used
;
; MODIFICATION HISTORY:
;	Written by:	Moritz Günther 21.4.2005
;	3.5.05		made ion_file and optional input
;-
pro TT_test_ioneq_rate, $
	ion_file=ion_file,$	;file name and FULL path with matrix with initial inonisation fractions of all ions in CHIANTi .ioneq format, one entry with t=t_ion
	outpath,	$	;Path and filename for output, without extension	
	abundfile=abundfile,$;elemental abundances, if not set cosmic.abund is used
	temperature=temperature,$
	result=result

common tt_units
common tt_constants
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
common tt_steps, delta_x, A, delta_t,x,stepcoeff
common tt_mean_atomic_weight,mean_atomic_weight

if n_params() lt 1 then begin
   printf, -2,"Necessary parameters missing" ;später ausführlicher machen
   return
endif

;--- Setup variables and calculate starting condition behind shock ---
print, "TT_SIMUL"
print, "Simulation of spectra from the accreting region of T Tauri stars"
print, ""
print, "Moritz Günther"
print, "Hamburger Sternwarte"
print, ""
a=1*cm^2
version='Test for tt_microscopic'
;initialize output
out=obj_new('tt_output',outpath,version)
out->log,'Start Run now'+systime()
out->log, 'Data will be stored at '+outpath
out->plog,'Test: Try to reproduce ioneq'
;set environmental variable 
if keyword_set (ion_file) then begin 
  read_ioneq, ion_file, ioneq_logt,ioneq,ioneq_ref 
  out->log,"Using Ionization equilib file:   "+ion_file
 endif else begin
  read_ioneq, concat_dir(!tt_data ,'test5.ioneq'),ioneq_logt,ioneq,ioneq_ref
  out->log,"Using Ionization equilib file:   "+concat_dir(!tt_data ,'test5.ioneq')
endelse  

if keyword_set(abundfile) then begin
  read_abund, abundfile ,abund,abund_ref
  out->log,"Using elemental abundances file: "+abundfile
 endif else begin
  read_abund, !abund_file ,abund,abund_ref
  out->log,"Using elemental abundances file: "+!abund_file
endelse  


TT_units
mean_atomic_weight=tt_get_mean_atomic_weight()

x_e=tt_get_xe(T_ion)
delta_x=1e+4		;set any non-zero starting value
delta_t=1e0*s

max_number_of_steps=3000

if(n_elements(temperature) eq 1) then begin 
  result=make_array(10,n_elements(ioneq[0,*,0]),n_elements(ioneq[0,0,*]))
  x_e_array=make_array(10)
  ;---hope to get equilib after many steps
  for number_of_steps=0,max_number_of_steps do begin
    if ((number_of_steps mod 100) eq 0) then print, 'temp:'+string(temperature)+'  step:'+string(number_of_steps)
    ;You can save a lot of time by disabeling bb_rad_loss
    tt_microscopic,temperature,temperature,1e+12,x_e,d_xe,en_loss,relabundHII,out, /rate_test
    ;check that equilibrium is constant
    if (max_number_of_steps-number_of_steps lt 10) then begin 
      result[max_number_of_steps-number_of_steps,*,*]=ioneq[0,*,*]
      x_e_array[max_number_of_steps-number_of_steps]=x_e
    endif  
  endfor
  if not(array_equal(result[0,*,*],result[1,*,*])) then begin
    ;Check, if difference is close to numerical precision (1e-7 for floats)
    index=where((abs(result[0,*,*]-result[1,*,*])/result[1,*,*] gt 1e-7) and result[0,*,*] gt 1e-4)
    if n_elements(index) gt 0 then begin
      temp=result[0,*,*]
      temp2=result[1,*,*]
      out->log,'After '+string(number_of_steps)+' no real equilib reached.'
      out->log,'absolute differences:'
      out->log,temp[index]-temp2[index]
    endif
  endif
  out->add_data, reform(ioneq[0,*,*]),temperature,[42]
  
endif else begin
  result=make_array(41,n_elements(ioneq[0,*,0]),n_elements(ioneq[0,0,*]))
  test=make_array(10,n_elements(ioneq[0,*,0]),n_elements(ioneq[0,0,*]))
  ;---loop over temperatures---
  for log_temp=40,80 do begin	;loop over all temprartures
    t_ion=10.^(float(log_temp)/10.)
    print,'Temperature:' ,t_ion,' K'
    ;---hope to get equilib after many steps
    for number_of_steps=0,max_number_of_steps do begin
      if ((number_of_steps mod 100) eq 0) then print, 'log_temp:'+string(alog10(t_ion))+'  step:'+string(number_of_steps)
      ;You can save a lot of time by disabeling bb_rad_loss in tt_microscopic
      tt_microscopic,t_ion,T_ion,1e+12,x_e,d_xe,en_loss,relabundHII,out, /rate_test
      ;check that equilibrium is constant
      if (max_number_of_steps-number_of_steps lt 10) then test[max_number_of_steps-number_of_steps,*,*]=ioneq[0,*,*]
    endfor
    if not(array_equal(test[0,*,*],test[1,*,*])) then begin
      ;Check, if difference is close to numerical precision (~1e-7 for floats)
      index=where((abs(test[0,*,*]-test[1,*,*])/test[1,*,*] gt 1e-7)and test[0,*,*] gt 1e-4)
      if (index[0] ne -1) then begin
        temp=test[0,*,*]
        temp2=test[1,*,*]
        out->log,'After '+string(number_of_steps)+' no real equilib reached.'
        out->log,'absolute differences:'
        out->log,temp[index]-temp2[index]
      endif
    endif
    out->add_data, reform(ioneq[0,*,*]),t_ion,[42,42,42,42,42]
    result[log_temp-40,*,*]=ioneq[0,*,*]
  endfor
  ;---plot results to check correctness----
  !p.multi=[0,2,2]
  ;---calculate ioneq with static condition
  print,'Calculate static solution'
  tt_test_ioneq, static
  ;compare to Mazzotta et al results
  read_ioneq,concat_dir(!xuvtop, 'ioneq/mazzotta_etal_ext.ioneq'),ioneq_logt,ioneq,ioneq_ref
  print,'dynamic (line) and static (dotted) results and Mazzotta et al (dash line) values are plotted'
  print,'if everything is fine, they are on top of each other'
  plot, (findgen(41)+40.)/10.,result[*,0,1],title='HII',Xtitle='log10(T in K)',yrange=[1e-5,1e0]
  oplot, (findgen(41)+40.)/10.,ioneq[*,0,1],line=2
  oplot, (findgen(41)+40.)/10.,static[*,0,1],line=1
  plot, (findgen(41)+40.)/10.,result[*,25,13],title='Fe XIV',Xtitle='log10(T in K)',/ylog,yrange=[1e-5,1e0]
  oplot, (findgen(41)+40.)/10.,ioneq[*,25,13],line=2
  oplot, (findgen(41)+40.)/10.,static[*,25,13],line=1
  plot,result[15,9,*],title='Ne at log(t)=5.5 and 7.0',xrange=[0,13],yrange=[1e-6,1e0],Xtitle='Ionisation stage',/ylog
  oplot,ioneq[15,9,*],line=2
  oplot,ioneq[30,9,*],line=2
  oplot,result[30,9,*]
  oplot,static[15,9,*],line=1
  oplot,static[30,9,*],line=1
  plot, (findgen(41)+40.)/10.,result[*,5,3],title='CIV',Xtitle='log10(T in K)',/ylog,yrange=[1e-5,1e0]
  oplot, (findgen(41)+40.)/10.,ioneq[*,5,3],line=2
  oplot, (findgen(41)+40.)/10.,static[*,5,3],line=1
endelse

out->write_data
;--- compare static solution and result of dynamic relaxation
;ignore where data is not given
test=abs(static[*,0:28,0:29]-result[*,0:28,0:29])
out->plog,'Maximum difference between static calculation and dynamic relaxation is'
out->plog,string(max(test))

end
;--- test results---
;These rate reproduce the Mazzotta et al. ioneq values
;Equilibrium is reached faster for heavier elements and higher temperatures.
;Only problems for H and He: needs far more than 300 1s steps to reach equilibrium
; -> 3000 steps -> equilibrium nearly reached
; one run takes about an hour this way, could be much less steps for higher temperatures
;but test in principle needs to run only once so it is a waste of working time to program a more 
;subtle routine
;for very long steps (like 1e4 s) possible problems with lsode solutions
;but this is far away from any useful application, so it was not investigated further



