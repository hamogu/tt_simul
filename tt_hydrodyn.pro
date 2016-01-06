;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_HYDRODYN
;
; PURPOSE:
;	This procedure solves the ODE for hydrodynamic movement of a plasma
;	The simulation adopts a two-fliud approch with an electronic component and
;	an ionic one. They move with the same bulk velocity but can have different
;	temperatures initially. The procedure simulates thhe energy flux between the components concentrating only on the
;	lightest ions, ionized hydogen. For a high degree of ionization collisions with neutral particales are unimportant and
;	less energy get transferred for heavier ions.
;	The equations are solved numerically with build-in IDL functions.
;
; CATEGORY:
;	Hydrodynamic
;
; CALLING SEQUENCE:
;	TT_HYDRODYN, t_ion,t_e,density,u,xe,en_loss,relativeabundHII
;
; INPUTS:
;	T_ion:		Ion temperature in Kelvin
;	t_e:		Electron temperature in Kelvin
;	density:	Density in g/cm^(-3)
;	u:		bulk velocity in cm/s 
;	xe:		number of free electrons per ion/atom
;	en_loss:	energy lost from the system (due to radiative processes/ionizations/...) 
;	relativeabundHII: ratio H-II/total number of heavy particles
;
; OUTPUTS:
;	T_ion:		Ion temperature in Kelvin
;	t_e:		Electron temperature in Kelvin
;	density:	Density in g/cm^(-3)
;	u:		bulk velocity in cm/s 
;
; KEYWORDS:
;	init:		if set, the conserved quantities (mass and momentum flux) are set,
;			no other calculation is done
;
; COMMON BLOCKS:
;	tt_Constants:	Needed to get Boltzmanns constant
;	tt_mean_atomic_weight: 	Well, obviously here for mean_atomic_weight
;	tt_steps:	for communication of step size delta_x and surface element A
;	tt_hydrodyn_internal: serves two purposes
;				1)keeps values of conserved quantities between calls
;				2)availability of variables throughout the functions of this module
;
; ADDITIONAL FUNCTIONS
;	This file consists of the main procedure and a few supplementary functions
;	Function taufunc,t_ion,t_e: calculates a generalized temperature
;	Function omega_eifunc, t_ion,t_e: calculates nergy trasfer from ions to electrons
;	Function velocity,t_ion,t_e: calculates velocity from temperature and momentum flux
;			Valid only for velocity < isothermal speed of sound
;	Function derivatives,x,y: calulates derivative of y=[t_ion,t_e] at x
;
; RESTRICTIONS
;	velocity < isothermal speed of sound
;	energy loss by radiation does not change over delta_x
;	xe does not change over delta_x
;
; ACCURACY
;	lsode should solve the ODEs to an absolute and relative error of 1e-07
;	assuming that the microscopic behaviour (radiation, xe) does not change
;
; TEST STATUS:
;	results give correct qualitative behaviour
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 15.03.2005
;	03.05.2005	modiefied to properly use the init keyword
;-

Function taufunc,t_ion,t_e
  common tt_hydrodyn_internal, j,p,mass,n_e,x_e,d_xe,Qcol,relabHII
  common tt_constants
  return, k_boltz*(t_ion+x_e*t_e)
end

Function omega_eifunc, t_ion,t_e
  common tt_hydrodyn_internal, j,p,mass,n_e,x_e,d_xe,Qcol,relabHII
  common tt_constants
  coulomblogarithm=9.4+1.5*alog(t_e)-0.5*alog(n_e)
  return, 3./2.*k_boltz*(t_ion-t_e)/t_e^(1.5)*coulomblogarithm/252.*relabHII
end
Function velocity,t_ion,t_e
  common tt_hydrodyn_internal, j,p,mass,n_e,x_e,d_xe,Qcol ,relabHII
  ;j^2 is too big for float, mass*j*j is OK
  return, p/(2.*mass*j)*(1.-sqrt(1.-4.*mass*j*j*taufunc(t_ion,t_e)/p^2.))
end

Function derivatives,x,y
  common tt_constants
  common tt_hydrodyn_internal, j,p,mass,n_e,x_e,d_xe,Qcol,relabHII
  t_ion=y[0]
  t_e=y[1]
  u=velocity(t_ion,t_e)
  tau=taufunc(t_ion,t_e)
  omega_ei=omega_eifunc(t_ion,t_e)
  k=k_boltz
  root=sqrt(p^2-4.*mass*j*j*tau)
  factor=-x_e*j/u
  numerator=omega_ei*(1./k+j*2.*tau/(3.*u*k*root))-t_ion*j*2.*Qcol/(root*3.*u)
  denominator=1.5*u+j*tau/root
  dt_ion= factor*numerator/denominator
  dt_e=t_e/t_ion*dt_ion+2.*j/(3.*k*u^2)*((1.+t_e/t_ion*x_e)*omega_ei-Qcol)-t_e*d_xe/x_e
  return, [dt_ion,dt_e]
end  


PRO TT_HYDRODYN, $
	t_ion,$
	t_e,$
	density,$	;density in g/cm^3
	u,$
	xe,$
	dxe,$
	en_loss,$
	relativeabundHII,$
	out,$
	Init=init

common tt_mean_atomic_weight
common tt_constants
common tt_steps, delta_x, A, delta_t,x

;local, stores conserved quantities between calls to ensure they are really conserved
common tt_hydrodyn_internal, j,p,mass,n_e,x_e,d_xe,Qcol,relabHII

n=density/mean_atomic_weight
n_e=xe*n

if (keyword_set(init) or n_elements(j) eq 0 or n_elements(p) eq 0) then begin
  j=n*u	;particle flux
  p=n_e*k_boltz*t_e+n*k_boltz*t_ion+n*u^2*mean_atomic_weight	;+x_e*n*u^2*m_e  ;momentum density flux
endif else begin 
  ;for testing
  out->log,'j:'+ string(j)+'   p:'+string(p)

  x_e=xe
  d_xe=dxe
  Qcol=en_loss
  relabHII=relativeabundHII
  mass=mean_atomic_weight	;+x_e*m_e ; m_ion=m_atom in mean_atomic_weight so in fact you should not add the electron mass

  ;check if velocity is subsonic
  ; sqrt(tau(t_ion,t_e)/mass) is isothermal speed of sound with "effectiv" mass and temperature
  if u gt sqrt(taufunc(t_ion,t_e)/mass) then out->plog,"ERROR in tt_hydrodyn - velocity is supersonic!"
  
  ;now take preferred ODE solver...
  y=[t_ion,t_e]
  status=1	;initialize every call - can be designed better, but for tests
  temp=lsode(y,x,delta_x,'derivatives',status)		;do integration x->x+delta_x
  if (status lt 1) then out->plog,'Numerival problems in integration, status os lsode:'+string(status)
  if (status lt -1) then out->plog,'ERROR in numerical integration'
  
  ;--- return values ---
  t_ion=temp[0]
  t_e=temp[1]
  u=velocity(t_ion,t_e)
  delta_t=delta_x/u		;u has changed!
  density=j/u*mean_atomic_weight
endelse
end
