;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_TEST_HYDRODYN
;
; PURPOSE:
;	This function run the tt_hydrodyn code in with different starting parameters and analyses
;	and displays the results. Printout informs about the expected analyical solution for theses test cases.
;
; CATEGORY:
;	Verification
;
; CALLING SEQUENCE:
;	TT_TEST_HYDRODYN,outpath
;
; INPUTS:
;	outpath:	path for output log file
;
;
; COMMON BLOCKS:
;	tt_units:	defines common unti symbols
;	tt_constants:	defines physical constants
;	tt_mean_atomic_weight: mean_atomic_weight in g
;	tt_steps:	step sizes
;
; EXAMPLE:
;		IDL> tt_test_hydrodyn,'./test_hyd'
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther, 2005-06-14
;	modified to new tt_hydrodyn keyword - no calculation if /init is set   Moritz Günther 14.07.2005
;-
pro TT_test_hydrodyn,outpath
common tt_units
common tt_constants
common tt_mean_atomic_weight, mean_atomic_weight
common tt_steps, delta_x, A, delta_t,x


;--- Setup variables and calculate starting condition behind shock ---
print, "TT_SIMUL"
print, "Simulation of spectra from the accreting region of T Tauri stars"
print, ""
print, "Moritz Günther"
print, "Hamburger Sternwarte"
print, ""

a=1*cm^2
version='Test for tt_hydrodyn'
;initialize output
out=obj_new('tt_output',outpath,version)
out->log,'Start Run now'+systime()
out->log, 'Data will be stored at '+outpath

TT_units
mean_atomic_weight=atomic_mass_unit

;setup counter, max_number_of_steps should never be reached but to avoid endless looping...
number_of_steps=0
max_number_of_steps=100 ;adjust later
x=0.
;---test steady flow---
out->plog,''
out->plog,'Testing steady flow'
;make data for test cases
u=		[1.   ,50.  ,100. ,120. ,10.  ,80.  ,80. ,80.  ,80.]*km/s
t_ion=		[3e6  ,3e6  ,3e6  ,3e6  ,3e5  ,3e7  ,3e6 ,3e6  ,3e6]
t_e=		[3e6  ,3e6  ,3e6  ,3e6  ,3e5  ,3e7  ,3e6 ,3e6  ,3e6]
density=	[1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-8,1e-16,1e-12]*g/cm^3
delta_x_test=	[1e-4 ,1e-4 ,1e-4 ,1e-4 ,1e-4 ,1e-4 ,1e-4,1e-4 ,1e-1]*s
test=[[t_ion],[t_e],[density],[delta_x_test]]
;in this simple case (T_ion=T_e) and no heat loss there should be steady flow
for j=0,n_elements(u)-1 do begin
  delta_x=delta_x_test[j]
  delta_t=delta_x/u[j]
  for i=0,10 do begin
    if (i eq 0) then tt_hydrodyn, t_ion[j],t_e[j],density[j],u[j],1.,0.,0.,1.,out,/init else tt_hydrodyn, t_ion[j],t_e[j],density[j],u[j],1.,0.,0.,1.,out
  endfor
endfor
if (max(abs(test-[[t_ion],[t_e],[density],[delta_x_test]])/test) lt 1e-5) then begin
  out->plog,'Steadflow OK within numerical precission'
endif else begin
  out->plog,'WARNING flow deviates from analytical solution!'  
  index=where(abs(test-[[t_ion],[t_e],[density],[delta_x_test]])/test)
  out->plog,'problems arise in situation number:'
  list=array_indices(test,index)
  out->plog, list[0,*]
  out->log,'index values of test array which differ from expectation'
  out->log,index
endelse  


;--- test supersonic flows ---
out->plog,' '
out->plog,'Testing supersonic flow warnings'
out->plog,'For 5 values there should be error messages for super sonic flow'
u=		[400. ,100. ,300. ,12.  ,100.]*km/s
t_ion=		[3e6  ,1e6  ,1e5  ,1e4  ,1e6]
t_e=		[3e6  ,0e0  ,1e7  ,0e4  ,0e6]
density=	[1e-12,1e-12,1e-12,1e-4 ,1e-12]*g/cm^3
delta_x_test=	[1e-4 ,1e-4 ,1e-4 ,1e-4 ,1e-2]*s
;there will be repeated numerical errors for supersonic flow, but there is no need to report them
;that is what the supersonic test is made for
oldexcept=!except
!except=0
for j=0,n_elements(u)-1 do begin
  delta_x=delta_x_test[j]
  delta_t=delta_x/u[j]
  tt_hydrodyn, t_ion[j],t_e[j],density[j],u[j],1.,0.,0.,1.,out,/init
  tt_hydrodyn, t_ion[j],t_e[j],density[j],u[j],1.,0.,0.,1.,out
endfor
out->plog,'End of supersonic flow test'
;--- test energy transfer from electron to ion gas
;keeping t_ion,u,density=const
t_ion=1e7
t_e=1e4
t_e_list=t_e
repeat begin
  delta_t=1e-1
  t_ion=1e7
  u=50.*km/s
  delta_x=u*delta_t
  density=1e-12*g/cm^3
  tt_hydrodyn, t_ion,t_e,density,u,1.,0.,0.,1.,out,/init ;keeping u,n const violates momentum conservation, so init with new values for every step
  tt_hydrodyn, t_ion,t_e,density,u,1.,0.,0.,1.,out
  t_e_list=[t_e_list,t_e]
endrep until (t_ion-t_e lt 5e5)
out->plog,'Keeping T_ion=const'
!p.multi=[0,2,2]  
plot, t_e_list, Xtitle='time in 1e-2 s',Ytitle='temperature' ,title='constant ion temperature'
 ;keeping t_e,u,density=const
 x=0
t_e=1e7
t_ion=1e4
t_ion_list=t_ion
x_list=0
repeat begin
  t_e=1e7
  u=50.*km/s
  density=1e-12*g/cm^3
  tt_hydrodyn, t_ion,t_e,density,u,1.,0.,0.,1.,out,/init
  tt_hydrodyn, t_ion,t_e,density,u,1.,0.,0.,1.,out
  t_ion_list=[t_ion_list,t_ion]
  x+=delta_x
  x_list=[x_list,x]
endrep until (t_e-t_ion lt 5e5)
out->plog,'Keeping T_electron=const'
plot, t_ion_list, Xtitle='time in  1e-2 s',Ytitle='temperature' ,title='constant electron temperature'
n=density/mean_atomic_weight
analytic=t_e*(1.-exp(-x_list*n*(9.4+1.5*alog(t_e)-.5*alog(n))/252./u/t_e^1.5))
oplot, analytic
out->plog, 'analytic solution superposed in plot'
;keeping u,density=const
t_ion=1e7
t_e=1e2
t_ion_list=t_ion
t_e_list=t_e
u=50.*km/s
density=1e-12*g/cm^3
u_list=u
repeat begin
  ;u=50.*km/s
  ;density=1e-12*g/cm^3
  tt_hydrodyn, t_ion,t_e,density,u,1.,0.,0.,1.,out,/init
  tt_hydrodyn, t_ion,t_e,density,u,1.,0.,0.,1.,out
  t_ion_list=[t_ion_list,t_ion]
  t_e_list=[t_e_list,t_e]
  u_list=[u_list,u]
endrep until (abs(t_e-t_ion) lt 5e5)
test=check_math()
!except=oldexcept
out->plog,'constant velocity and density'
out->plog,'lines should be symmetric to 5e+06'
if max(abs((t_ion_list+t_e_list)-1e7)) lt 1e3 then out->plog,'Symmetrical within numerical precision' 
plot, t_ion_list, Xtitle='time in 1e-2 s',Ytitle='temperature' ,title='constant velocity and density',Yrange=[1e4,1e7]
oplot, t_e_list
end
