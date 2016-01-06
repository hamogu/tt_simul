;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_UNITS
;
; PURPOSE:
;	This procedure set units and constants of nature which are avaiable
;	to other routines by referencing the common block needed
;
; CATEGORY:
;	mixed
;
; CALLING SEQUENCE:
;	TT_UNITS
;
; DEFINITION COMMON BLOCK:
;	TT_Units:	units are transfered to the cgs system, e.g. 1m=100cm
;	TT_Constats:	phsical constants
;			electron mass, proton mass, k_boltz, u, mass_number
;			for elemts up to Zn
; BIBLOGRAPHY:
;	Unsöld/Baschek: "Der neue Kosmos" 5. Auflage
;	www.webelements.com
;	Kuchling, "Taschenbuch der Physik"
;
;
; TEST STATUS:
;	checked
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther 20.01.2005
;-


pro TT_UNITS
common TT_units, km,m,cm,mm,s,minute,kg,g,erg,kelvin,ly,pc,ae,angstroem,year
common TT_constants, m_e,m_proton,m_H, k_boltz,mass_number,atomic_mass_unit,e_charge,c_light,h_planck,r_sun,m_sun
cm=1.
m=100.*cm
km=1000.*m
mm=cm/10.
s=1.
minute=60.*s
g=1.
kg=1000.*g
erg=1.
Kelvin=1.
AE=1.496e11*m
au=ae
year=3.1558e7*s
ly=9.4606e12*km
pc=3.086e16*m  
angstroem=1e-10*m  
m_e=9.1094e-28*g
m_proton=1.6726e-24*g
m_H=1.6736e-24*g
k_Boltz=1.38066e-16*erg/Kelvin
atomic_mass_unit=1.66054e-24*g
mass_number=[1.0,4.0,6.94,9.01,10.81,12.01,14.01,16.0,19.0,20.18, $
	22.99,24.31,26.98,28.09,30.93,32.07,35.45,39.95,$
	39.10,40.08,44.96,47.87,50.94,52.0,54.04,55.85,58.93,58.69,63.55,65.41]
e_charge=1.60217646e-19 ;Coulomb
c_light=2.9979e+08*m/s
h_planck=6.6261e-27*erg*s
r_sun=6.96e8*m
m_sun=1.989d33*g
end
