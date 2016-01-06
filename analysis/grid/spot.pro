;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	spot
;
; PURPOSE:
;	This procedure calculates spot charcteristics: size, filling factor and mass flux.
;
; CATEGORY:
;	analyses
;
; CALLING SEQUENCE:
;	 spot,lines, index_min, obsfile,r_star=r_star, d_star=d_star, simdata_file=simdata_file,abund_factor=abund_factor
;
; INPUTS:
;	lines:		a string array which specifies which lines should be used
;	index_min:	1-dim index of the best-fit model
;	obsfile:	path to observational data in a format compatible to read_obs_data.pro

; OPTIONAL INPUTS:
;	r_star:		radius of the star relative to the solar radius, default=1
;	d_star:		distance to the star in pc, default=57pc (TW HYA)
;	simdata_file:	path to an IDL .sav with with a data structure "data" which contains an output of analyse_grid.pro
;			default is rad6e3newfe.sav in the working directory
;	abund_factor:	array with elemental abundances, size= [number fo lines,2]
;			abund_factor[*,0] is the factor compared to Allen(1975), e.g. from fit_abund
;			abund_factor[*,1] are the corresponding errors
;	nh		hydrogen column density in cm^-2 which is used for correction of the observations
;
; OUTPUT:
;	A_spot, filling factor and mass flux on screen
;
; KEYWORDS:
;	ratiofile:	if set ratios is not a string array, but the path to a file with such an array written line by line	
;
; TEST STATUS:
;	tested
;
; MODIFICATION HISTORY:
;	Written by:	Moritz Guenther 23.09.2005
;-


pro spot,lines, index_min, obsfile,r_star=r_star, d_star=d_star, simdata_file=simdata_file,abund_factor=abund_factor,nh=nh

common tt_units
common tt_constants

if keyword_set(simdata_file) then begin
  restore,simdata_file
  simdata=data
endif else begin 
  restore,'rad6e3newfe.sav'
  simdata=data
endelse  

n_lines=n_elements(lines)

r_star= n_elements(r_star) eq 0 ? 6.96e8*m : r_star*6.96e8*m  ;relativ to our sun's size
d_star= n_elements(d_star) eq 0 ? 57d*pc : d_star*pc*1d0      ;convert to double
 
read_obs_data, obsfile,obsdata,obs_ref,nh=nh

a_spot=make_array(n_lines,2,/double)

;if not given, set abund_factor to 1 for each line and all errors ro 0
if n_elements(abund_factor) eq 0 then abund_factor=[[1],[0]]  ;this line in necessary to make the next line legal, if abund_factor is formerly undefined
if n_elements(abund_factor[*,0]) ne n_lines then begin
  if n_elements(abund_factor) ne 0 then print, 'abund_factor does not match lines ->setting to 1'
  abund_factor=make_array(n_lines,2,val=1.)
  abund_factor[*,1]=0.	;no errors
endif  

for i=0,n_lines-1 do begin
  index=where(strlowcase(lines[i]) eq strlowcase(obsdata.name))
  ;lines[i] should be the name of a line in the simdata structure
  test=execute('sim=(simdata.'+lines[i]+')[index_min]')
  if test ne 1 then begin
    print,'Line:',lines[i],' does not exist or index_min is out of range'
    print,'aborting'
    return
  endif  
  
  a_spot[i,0]=(obsdata.flux)[index]/sim*4.*!pi*d_star^2./abund_factor[i,0]
  a_spot[i,1]=1./sim*4.*!pi*d_star^2.*sqrt(((obsdata.error)[index]/abund_factor[i,0])^2.+((obsdata.flux)[index]/abund_factor[i,0]^2.*abund_factor[i,1])^2)
  
endfor

f=a_spot/(4*!pi*r_star^2.) ;filling factor

m_sun=1.989e30*kg
year=3.1157e7*s
;mass infall in solar masses per year
M_dot=a_spot*(simdata.n0)[index_min]*m_h*1.245/m_sun *(simdata.v0)[index_min]* year  ;1.245 is mean molecular weight

;weighted mean
a_spot_mean=total(a_spot[*,0]/a_spot[*,1]^2.)/total(1./a_spot[*,1]^2.)
;2 ways to get the error
a_spot_sigma_gauss=sqrt(1./total(1./(a_spot[*,1]^2.))) ;error propagation by usual laws
a_spot_sigma_direct=n_lines le 1 ?0:sqrt(total((a_spot[*,0]-a_spot)^2./a_spot[*,1]^2.)/(float(n_lines-1.)*total(1./a_spot[*,1]^2.)))  ;error as seen from the basic sample
a_spot_sigma=max([a_spot_sigma_gauss,a_spot_sigma_direct])  ;take whcichever is higher

f_mean=a_spot_mean/(4*!pi*r_star^2.)
f_mean_sigma=a_spot_sigma/(4*!pi*r_star^2.)

M_dot_mean=a_spot_mean*(simdata.n0)[index_min]*m_h*1.245/m_sun *(simdata.v0)[index_min]* year
M_dot_sigma=a_spot_sigma*(simdata.n0)[index_min]*m_h*1.245/m_sun *(simdata.v0)[index_min]* year

print, 'Calculate spot characteristics'
print, 'Stellar parameters'
print, 'R_star=',r_star/(6.96e8*m),' R_sun'
print, 'distance=',d_star/pc,' pc'
print,format='("v0=",I3," [km/s],  log(n0)=",f6.1," log([1/cm^3])")',simdata.v0[index_min]/km,alog10(simdata.n0[index_min])
print, '--------------------------------------------------------------------------------------------------------------------'
print, '            line     A_spot[cm^2]          sigma               f           sigma       M_dot[M_sun/year]       sigma'
															
for i=0,n_lines-1 do print, format='((a16)6(e16.2))',lines[i],a_spot[i,0],a_spot[i,1],f[i,0],f[i,1],m_dot[i,0],m_dot[i,1]
print, '--------------------------------------------------------------------------------------------------------------------'
print, format='("MEAN RESULTS    "6(e16.2))',a_spot_mean,a_spot_sigma,f_mean,f_mean_sigma,m_dot_mean,m_dot_sigma

end