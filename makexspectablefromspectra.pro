;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	MAKEXSPECTABLEFROMSPECTRA
;
; PURPOSE:
;	This procedure produces a *.fits file to be read into XSPEC as a table model.
;	It reads sperctra created from a rectangular grid of TT_simul runs. The spectral files should be called
;	e.g. 3_5spectrum.sav and the values n0 and vo need to adhere to the TT_Simul naming convention.
;
; CATEGORY:
;	Analysis
;
; CALLING SEQUENCE:
;	MAKEXSPECTABLEFROMSPECTRA,intpar1,intpar2,naddparm,fitsfile,modelname,addnames,d=d,r=r
;
; INPUTS:
;	intpar1:	Array of numbers specifying the v0, e.g. [0,2,4] for 0_xspectrum.sav...
;	intpar1:	Array of numbers specifying the n0, e.g. [0,2,4] for x_0spectrum.sav...
;	naddparam:	Number of additional paramters, here elements with one spectra, typical naddparam=8
;	fitsfile:	output file name
;	modelname:	Model name shown in XSPEC
;	addnames:	name of additional params in XSPEC, typical ; addnames=['C','N','O','Ne','Mg','Si','S','Fe']
;
; OPTIONAL INPUTS
;	d:		distance to star on pc [default=10pc]
;	r:		radius of star in solar radii
;		If d and r are given then the norm in XSPEC is later equvalent to the surface filling factor.
;		If the r i not given, then the norm will show the A_Spot in untis of 10^20 cm^2.
;
; EXAMPLE
;	makexspectablefromspectra,[8,9,10,11,12,13,14,15],[2,3,4,5,6,7,8],8,'name.fits','V4046 Sgr', ['C','N','O','Ne','Mg','Si','S','Fe'], 83, 1.2
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 2006/2007
;	29.08.2008	Modified the normalisation, added HDUCLASS header keywords
;-




pro makexspectablefromspectra,intpar1,intpar2,naddparm,fitsfile,modelname,addnames,d=d,r=r

if n_elements(d) ne 1 then d=10. ;default 10 pc distance


common tt_units



;--- generate wavelength and energy scales ---
min_lambda=1.
max_lambda=50.
ang2kev=12.39854
kev_scale=findgen(12000)*0.00025+0.000125
lambda_scale=ang2kev/kev_scale
lambda_scale=lambda_scale[where(lambda_scale ge min_lambda and lambda_scale le max_lambda)]
n_energies=11008

binboundaries=mid2bound(lambda_scale)
upper=binboundaries[0:n_elements(lambda_scale)-1]
lower=binboundaries[1:*]
bin_width=upper-lower
e_max=ang2kev/lower	;lower wavelength->higher energy!
n_energies=n_elements(e_max)
e_min=[0,e_max[0:n_energies-2]]

;--- Setup parameters

nintparm=2
intpoints=n_elements(intpar1)*n_elements(intpar2)

v0=(findgen(40)*25.+200.)
n0_index=indgen(22)-9
n0=10.^(findgen(22)/3.+7.);* 2.08316e-24


;v0=(findgen(17)*25.+200.)
;n0=10.^(findgen(13)/3.+10.)
v0=v0[intpar1]
n0=n0[intpar2+9]

;--- Fits header ---
fxhmake,header,/initi,/extend,/date
fxaddpar,header, 'MODLUNIT','photons/cm^2/s','Unit of spectrum'
fxaddpar,header,'MODLNAME',modelname
fxaddpar,header,'ORIGIN','TT_Simul','by Moritz Guenther, Hamburger Sternwarte'
fxaddpar,header,'REDSHIFT','F','redshift is not a fit parameter'
fxaddpar,header,'ADDMODEL','T','additive model'
fxaddpar,header,'HDUCLASS', 'OGIP    '
fxaddpar,header,'HDUDOC  ', 'OGIP/92-009'
fxaddpar,header,'HDUCLAS1', 'XSPEC TABLE MODEL'
fxaddpar,header,'HDUVERS1', '1.0.0   '

fxaddpar,header, 'D_star',d,'distance to object used to scale fluxis in pc'
if n_elements(r) eq 1 then begin
  fxaddpar,header,'R_star','radius of star in R_sun'
  fxaddpar,header,'unitnorm','norm in XSPEC represents the filling factor'
endif else begin
  fxaddpar,header,'unitnorm','norm in XSPEC represents A_spot [10^20 cm^2]'
endelse


;--- Parameters section ---
fxaddpar,header1,'EXTNAME','PARAMETERS'
fxaddpar,header1,'NINTPARM',2
fxaddpar,header1,'NADDPARM',naddparm
fxaddpar,header1,'HDUCLASS', 'OGIP    '
fxaddpar,header1,'HDUCLAS1', 'XSPEC TABLE MODEL'
fxaddpar,header1,'HDUCLAS2', 'PARAMETERS'
fxaddpar,header1,'HDUVERS1', '1.0.0   '


parameters=replicate({NAME:'v0          ',METHOD:long(0),INITIAL:525.,DELTA:25.,MINIMUM:min(v0),BOTTOM:min(v0),TOP:max(v0),MAXIMUM:max(v0),NUMBVALS:long(n_elements(v0)),VALUE:v0},naddparm+nintparm)
;n0
parameters[1].name='n0'
parameters[1].initial=1e12
parameters[1].delta=1e11
parameters[1].minimum=min(n0)
parameters[1].bottom=min(n0)
parameters[1].top=max(n0)
parameters[1].maximum=max(n0)
parameters[1].numbvals=long(n_elements(n0))
parameters[1].value=n0
; additional parameters
for i=nintparm,naddparm+nintparm-1 do begin
  parameters[i].name=addnames[i-nintparm]
  parameters[i].initial=1.
  parameters[i].delta=1.
  parameters[i].minimum=0.
  parameters[i].bottom=.1
  parameters[i].top=3.
  parameters[i].maximum=10.
  parameters[i].numbvals=long(1)
endfor

;--- Energy section ---
fxaddpar,header2,'EXTNAME','ENERGIES'
fxaddpar,header2,'HDUCLAS1', 'XSPEC TABLE MODEL'
fxaddpar,header2,'HDUCLAS2', 'ENERGIES'
fxaddpar,header2,'HDUVERS1', '1.0.0   '


energies=replicate({ENERG_LO:0.,ENERG_HI:0.},n_energies)


;--- Spectra section ---
fxaddpar,header3,'EXTNAME','SPECTRA'
fxaddpar,header3,'HDUCLAS1', 'XSPEC TABLE MODEL'
fxaddpar,header3,'HDUCLAS2', 'MODEL SPECTRA'
fxaddpar,header3,'HDUVERS1', '1.0.0   '


basespectra={PARAMVAL:[0.,0.],INTPSPEC:fltarr(n_energies)}
for i=1,naddparm do basespectra=add_tag(basespectra,fltarr(n_energies),'ADDSP00'+strtrim(string(i),2))
;spectra=replicate({PARAMVAL:0.,INTPSPEC:fltarr(n_energies),ADDSP001:fltarr(2350),ADDSP002:fltarr(2350),ADDSP003:fltarr(2350),ADDSP004:fltarr(2350),ADDSP005:fltarr(2350),ADDSP006:fltarr(2350),ADDSP007:fltarr(2350),ADDSP008:fltarr(2350)},intpoints)
spectra=replicate(basespectra,intpoints)


for ipar1=0,n_elements(intpar1)-1 do begin  ;cannot use i here, because that is also contained in the .sav files which makes trouble
  for jpar2=0,n_elements(intpar2)-1 do begin
     restore,strcompress(string(intpar1[ipar1])+'_'+string(intpar2[jpar2]),/remove_all)+'spectrum.sav'
     spectra[ipar1*n_elements(intpar2)+jpar2].paramval= base_spectra.paramval/[km,1.]
  endfor
endfor

factor=1.

;CAREFUL: Definition of units changed 29.08.2008!

if n_elements(r) eq 1 then begin
  r=r*6.961e10  ;solar radii
  factor=r*r/d/d/pc/pc
endif else begin
  factor=1e20/(4.*!pi*d*pc)/(d*pc) ; need to break up pc*pc to avoid floating underflow 
endelse

; The spectra are given in photons/s/Ang/sr/cm^2, but XSPEC expects photons/s/cm^2
 factor=factor*4.*!pi*bin_width


for ipar1=0,n_elements(intpar1)-1 do begin
  for jpar2=0,n_elements(intpar2)-1 do begin
     restore,strcompress(string(intpar1[ipar1])+'_'+string(intpar2[jpar2]),/remove_all)+'spectrum.sav'
     ; in some simulations there are nan values in the spectrum
     ; interpolate over those regions, but give warning
     for k=0,naddparm do begin
         indnan = where(base_spectra.(k+1) ne base_spectra.(k+1))
         ind = where(base_spectra.(k+1) eq base_spectra.(k+1))
         if indnan[0] ne -1 then begin
             print, 'Interpolate over '+string(n_elements(indnan))+' nan values in spectrum'+strcompress(string(intpar1[ipar1])+'_'+string(intpar2[jpar2]),/remove_all)
             n_elem = findgen(n_elements(base_spectra.(k+1)))  ; simply use element number as x axis, could use lambda or energy, but this should do
             base_spectra.(k+1)[indnan] = spline(n_elem[ind], base_spectra.(k+1)[ind],n_elem[indnan]) 
         endif
     endfor
     spectra[ipar1*n_elements(intpar2)+jpar2].intpspec= base_spectra.intpspec*factor
     for k=0,naddparm-1 do spectra[ipar1*n_elements(intpar2)+jpar2].(k+2)=base_spectra.(k+2)*factor
  endfor
endfor


for i=0,n_energies-1 do energies[i].energ_lo=e_min[i]
for i=0,n_energies-1 do energies[i].energ_hi=e_max[i]


mwrfits,test,fitsfile,header,/create
mwrfits,parameters,fitsfile,header1
mwrfits,energies,fitsfile,header2
mwrfits,spectra,fitsfile,header3

end
