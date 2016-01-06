;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_SIMUL
;
; PURPOSE:
;	This procedure is the main program of the TT_SIMUL project.
;	It initialises the output, and orginises the physical simulation, firt the shock accoding
;	to the Rankine-Hugoniot conditions  then step wise the post shock zone.
;
; CATEGORY:
;	MAIN
;
; CALLING SEQUENCE:
;	TT_SIMUL,u,density,T_ion,T_e,ionfrac_file,outpath,T_star=T_star,abundfile=abundfile,output=output,paramfile=paramfile
;
; INPUTS:
;	u:		bulk velocity in cm/s 
;	density:	Density in g/cm^(-3)
;	T_ion:		Ion temperature in Kelvin
;	t_e:		Electron temperature in Kelvin
;	ionfrac_file:	*.ioneq file according to CHIANTI specifiacations with ioneq values
;			for a single temperature as pre-shock starting point
;	outpath:	path for output files
;
; OPTIONAL INPUTS:
;	T_star:		temperature of a blackbody background field for level population calculations, default is 6000K
;	abundfile:	CHIANTI *.abund file, default is !abund
;	paramfile:	path to a file with program parameters, default is './tt_simul.par'
;
;  KEYWORDS:
;	quiet		suppresses status messages after every step
;
; OPTIONAL OUTPUTS:
;	output:		reference to output object for further direct analyses
;
; COMMON BLOCKS:
;	tt_units:	unit signs (km, m,cm,...)
;	tt_Constants:	
;	tt_mean_atomic_weight: 	Well, obviously here for mean_atomic_weight
;	tt_steps:	for communication of step size delta_x and surface element A
;	tt_elements:	ioneq and abund from CHIANTI
;
; TEST STATUS:
;	results give correct qualitative behaviour
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 15.05.2005
;	20.06.2005	Moritz Guenther - fixed memory leak (out was not destroyed)
;	14.07.2005	Moritz Guenther - added quiet keyword
;	17.10.2006	Moritz Guenther - remove the execute command by including standard values for all input_params and add input keywords
;-
pro TT_SIMUL, u, $
	density,$ 	;density in g/cm^-3 of infalling material
	T_ion,$		;Temperature in K of infalling ions
	T_e,$		;Temperature in K of infalling electrons
	ion_frac_file,$	;file name and FULL path with matrix with initial inonisation fractions of all ions in CHIANTi .ioneq format, one entry with t=t_ion
	outpath,$	;Path and filename for output, without extension
	T_star=T_star,$	;Temperature of background radiation, calculation stops when gas is cooled down to T_star		
	abundfile=abundfile,$;elemental abundances, if not set cosmic.abund is used
	output=output,$	;returns the output structure
	$;paramfile=paramfile,$; file with run parameters
	quiet=quiet,	$;suppresses some printout
	max_number_of_steps=max_number_of_steps,$  	;number of steps before program terminates (e.G. in case of endless loops)
	max_depth=max_depth,$				;depth were program cancels calculation (e.g. in case of unsuitable input parameters)
	max_hydrodyn_change=max_hydrodyn_change,$	;max. relative change in hydrodynamic variables per step
	max_ioneq_change=max_ioneq_change,$		;max. relative change in any ionisation state per step
	min_stepsize=min_stepsize			;minimum stepsize in case of steep gradients

common tt_units
common tt_constants
common tt_mean_atomic_weight, mean_atomic_weight
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
common tt_steps, delta_x, A, delta_t,x

starttime=systime(/seconds)

version='1.1'

if n_params() lt 6 then begin
   printf, -2,"Necessary parameters missing" 
   return
endif

;initialize output
out=obj_new('tt_output',outpath,version)
out->log,'Start Run now'+systime()
out->log, 'Data will be stored at '+outpath
output=out
; 
; ;read parameter file
; if n_elements(paramfile) eq 0 then paramfile='tt_simul.par'
; openr,lu,paramfile,/get_lun
; command='' ;to define it as a string
; while (not eof(lu)) do begin
;   readf, lu, command
;   ;print, command
;   if not execute(command) then print,'ERROR'
;   out->log,command
; endwhile
; free_lun, lu






;--- Setup variables and calculate starting condition behind shock ---
print, "TT_SIMUL"
print, "Simulation of spectra from the accreting region of T Tauri stars"
print, ""
print, "Moritz Guenther"
print, "Hamburger Sternwarte"
print, ""

;set environmental variable 
defsysv,'!ioneq_file', ion_frac_file
out->log,"Using Ionization equilib file:   "+!ioneq_file
if keyword_set(abundfile) then defsysv, '!abund_file', abundfile
out->log,"Using elemental abundances file: "+!abund_file

read_abund, !abund_file ,abund,abund_ref
read_ioneq, !ioneq_file ,ioneq_logt,ioneq,ioneq_ref

TT_units
mean_atomic_weight=tt_get_mean_atomic_weight()

;setup counter, max_number_of_steps should never be reached but to avoid endless looping...
number_of_steps=0

if not keyword_set(T_star) then T_star=6000.
x_e=tt_get_xe(T_ion)
x=0.
delta_x=1e5		;set any non-zero starting value
delta_t=delta_x/u
a=1*cm^2
if not keyword_set(max_number_of_steps) then max_number_of_steps=10000 
if not keyword_set(max_depth) then max_depth=1e5*km
if not keyword_set(max_hydrodyn_change) then max_hydrodyn_change=0.05
if not keyword_set(max_ioneq_change) then max_ioneq_change=0.05
if not keyword_set(min_stepsize) then min_stepsize=1e2*cm



;Temperature in ioneq is very high 1e+42 to ensure that this value will not be used for later interpolation
out->add_data, reform(ioneq[0,*,*]),1e+42,[-1,density,density/mean_atomic_weight,u,t_ion,t_e,x_e,0.]

;in the shock (approximatied by a mathematical discontinuity) all veriables change
;according to the Rankine-Hugoniot jump conditions
; The electrons are only affected by external forces and just get compressed adiabatically
old_density=density
TT_RANKINE_HUGONIOT_TURHO, t_ion,u,density		;for ions
T_e=T_e*(density/old_density)^(5./3.-1.)		;adiabatic compression
delta_t=delta_x/u					;u has changed!
number_of_steps++
out->add_data, reform(ioneq[0,*,*]),1e+41,[0,density,density/mean_atomic_weight,u,t_ion,t_e,x_e,0.]

tt_hydrodyn, t_ion,t_e,density,u,x_e,d_xe,en_loss,relabundHII,out,/init
;---Main loop--
repeat begin
  if keyword_set(quiet) then begin
    out->log,"Calculating step "+string(number_of_steps)+" in depth "+string(x/km)+" km"
    out->log,'Try delta_x='+string(delta_x/1e5)+' km'
  endif else begin
    out->plog,"Calculating step "+string(number_of_steps)+" in depth "+string(x/km)+" km"
    out->plog,'Try delta_x='+string(delta_x/1e5)+' km' 
  endelse  
  ;--- keep old values ---
 ;on first step there are no old values for d_xe and variables right of it
  old_hydrodyn=n_elements(d_xe) eq 0 ? [density,u,t_ion,t_e,x_e]:[density,u,t_ion,t_e,x_e,en_loss,relabundHII]
  old_ioneq=ioneq
  ;---try step
  tt_microscopic, t_ion,T_e,density/mean_atomic_weight,x_e,d_xe,en_loss,relabundHII,out,radtemp=t_star
  tt_hydrodyn, t_ion,t_e,density,u,x_e,d_xe,en_loss,relabundHII,out
  ;---test if step is too large ---
  ;on first step there are no old values for d_xe and variables right of it
  new_hydrodyn=n_elements(d_xe) eq 0 ? [density,u,t_ion,t_e,x_e]:[density,u,t_ion,t_e,x_e,en_loss,relabundHII]
  hydrodyn_change= max(abs(old_hydrodyn-new_hydrodyn)/new_hydrodyn)
  ioneq_change=max(abs(old_ioneq-ioneq))
  if((hydrodyn_change gt max_hydrodyn_change or ioneq_change gt max_ioneq_change)and delta_x gt min_stepsize) then begin
    ;reset and redo step with smaller stepsize
    ioneq=old_ioneq
    density=old_hydrodyn[0]
    u=old_hydrodyn[1]
    t_ion=old_hydrodyn[2]
    t_e=old_hydrodyn[3]
    x_e=old_hydrodyn[4]
    ;the rest are output parameters of tt_microscopic, there is no need to reset them
    endif else begin
    ;keep step
     x+=delta_x
    out->add_data, reform(ioneq[0,*,*]),t_e*(1+t_ion/t_e*m_e/mean_atomic_weight),[x,density,density/mean_atomic_weight,u,t_ion,t_e,x_e,en_loss]
    ;out->write_data
    number_of_steps++
  endelse
  ;---adjust step size---  
  scale=hydrodyn_change/(0.5*max_hydrodyn_change) > ioneq_change/(0.5*max_ioneq_change)
  delta_x/=scale
  delta_t/=scale
  ;fluid approximation is not valid for small steps
  if delta_x lt min_stepsize then begin
    delta_x=min_stepsize
    delta_t=delta_x/u
  endif  
endrep until (T_ion le (2*T_star>1e4) || T_e le (2*T_star>1e4) || number_of_steps ge max_number_of_steps || x gt max_depth ||$
 ~ array_equal([x,density,u,t_ion,t_e,x_e],[x,density,u,t_ion,t_e,x_e]) ) ;checking for NaN's; The Chianti functions are meant for temperatures > 10.000 K, so wahtever the stellar temperature is, we need to stop at latest there.
;--- Finishing ---
out->write_data
out->log, "Program stoped after "+string(number_of_steps)+" intervals"
out->log, "in depth "+string(x/km)+" km"
if number_of_steps ge max_number_of_steps then out->log, "Programm aborted because max. number of steps reached"
if ~ array_equal([x,density,u,t_ion,t_e,x_e],[x,density,u,t_ion,t_e,x_e]) then out->plog,'ERROR - NaN detected - aborting'

out->plog,'Simulation took '+string(systime(/seconds)-starttime)+' sec'

; cleanup and free memory
if n_elements(output) eq 0 then begin
  obj_destroy,out
endif  

end
