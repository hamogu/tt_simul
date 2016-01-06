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
;	en_loss:	energy lost from the system (due to radiative processes/ionizations/...) in später
;	relativeabundHII: ratio H-II/total number of heavy particles
;
; OUTPUTS:
;	T_ion:		Ion temperature in Kelvin
;	t_e:		Electron temperature in Kelvin
;	density:	Density in g/cm^(-3)
;	u:		bulk velocity in cm/s 
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
; 	Written by:	Moritz Günther 15.03.2005
;-

Function tt_hydrodyn::taufunc,t_ion,t_e
  common tt_constants
  return,k_boltz*(t_ion+self.x_e*t_e)
end

Function tt_hydrodyn::omega_eifunc, t_ion,t_e
  common tt_constants
  coulomblogarithm=9.4+1.5*alog(t_e)-0.5*alog(self.x_e*self.n)
  return, 3./2.*k_boltz*(t_ion-t_e)/t_e^(1.5)*coulomblogarithm/252.*self.relabHII
end
Function tt_hydrodyn::velocity,t_ion,t_e
  ;j^2 is too big for float, m*j*j is OK
  return, self.p/(2.*self.m*self.j)*(1.-sqrt(1.-4.*self.m*self.j*self.j*self->taufunc(t_ion,t_e)/self.p^2.))
end

Function tt_hydrodyn::derivatives,x,y
  common tt_constants
  common tt_hydrodyn_internal, j,p,m,n_e,x_e,Qcol,relabHII
  t_ion=y[0]
  t_e=y[1]
  u=self->velocity(t_ion,t_e)
  tau=self->taufunc(t_ion,t_e)
  omega_ei=omega_eifunc(t_ion,t_e)
  k=k_boltz
  root=sqrt(self.p^2-4.*self.m*self.j*self.j*tau)
  factor=-self.x_e*self.j/self.u
  numerator=omega_ei*(1./k+self.j*2.*tau/(3.*u*k*root)-t_ion*self.j*2.*self.en_loss/(root*3.*u))
  denominator=1.5*u+self.j*tau/root
  dt_ion= factor*numerator/denominator
  dt_e=t_e/t_ion*dt_ion+2.*self.j/(3.*k*u^2)*((1.+t_e/t_ion*self.x_e)*omega_ei-self.en_loss)
  return, [dt_ion,dt_e]
end  


PRO TT_HYDRODYN::calculate_step,x_e,en_loss,relabHII,delta_x
  common tt_constants
self.x_e=x_e
self.en_loss=en_loss
self.relabHII=relabHII

;check if velocity is subsonic
; sqrt(tau(t_ion,t_e)/m) is isothermal speed of sound with "effectiv" mass and temperature
if self.u gt sqrt(self->taufunc(self.t_ion,self.t_e)/self.m) then out->plog,"ERROR in tt_hydrodyn - velocity is supersonic!"

;now take preferred ODE solver...
y=[self.t_ion,self.t_e]
temp=lsode(y,self.x,delta_x,'tt_hydrodyn::derivatives',self.status)		;do integration 0->delta_x
						;but could do x->x+delta_x as well
;print, "Status of ODE integration: "+string(status)
if self.status ne 2 then self.out->plog,'Status of lsode after calculating is'+string(self.status)
;--- return values ---
self.t_ion=temp[0]
self.t_e=temp[1]
self.u=velocity(self.t_ion,self.t_e)
self.x+=delta_x
end

function tt_hydrodyn::INIT,mean_atomic_weight,n,u,x_e,t_e,t_ion,pout,x
  common tt_constants
  self.m=mean_atomic_weight
  self.j=n*u
  self.p=n*x_e*k_boltz*t_e+n*k_boltz*t_ion+n*u^2*mean_atomic_weight+x_e*n*u^2*m_e
  self.out=pout
  self.t_ion=t_ion
  self.t_e=t_e
  self.u=u
  self.n=n
  self.x=x
  return,1
end
pro  tt_hydrodyn::CLEANUP
  IF (PTR_VALID(self.out)) THEN PTR_FREE,self.out
end
function tt_hydrodyn::u
  u=self.u
  return, u
end
function tt_hydrodyn::t_ion
  t_ion=self.t_ion
  return, t_ion
end
function tt_hydrodyn::t_e
  t_e=self.t_e
  return, t_e
end    
function tt_hydrodyn::n
  n=self.n
  return, n
end 
function tt_hydrodyn::rho
  rho=self.n*self.m
  return, rho
end
 function tt_hydrodyn::x
  x=self.x
  return, x
end
pro tt_hydrodyn__define
  void={tt_hydrodyn,j:0.,p:0.,m:0.,out:obj_new(),u:0.,t_ion:0.,t_e:0.,n:0.,status:1,x:0.,x_e:0.,relabHII:0.,en_loss:0.}
end
 
