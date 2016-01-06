;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	fit
;
; PURPOSE:
;	This procedure performs a Chi-square fit from observational data on a grid of models. The parameters for the best-fit model are printed 
;	and a contourplot with the chi-squared values. Because the model depends non-linearly on the parameters n0 and v0 and the meassurement errors
;	are non-gaussian (They are poisson distributed) no simple analytical formula gives the relation of probabilities to chi squared.
;
; CATEGORY:
;	fitting
;
; CALLING SEQUENCE:
;	FIT, ratios, datafile,lev=lev,ratiofile=ratiofile
;
; INPUTS:
;	ratios:		a string array with the ratios to be computed in a format compatible to calculate_ratios.pro
;	datafile:	path to observational data in a format compatible to read_obs_data.pro
;
; OPTIONAL INPUTS:
;	simdata_file:	path to an IDL .sav with with a data structure "data" which contains an output of analyse_grid.pro
;			default is rad6e3newfe.sav
;	nh		hydrogen column density in cm^-2 which is used for correction of the observations
;
; KEYWORDS:
;	ratiofile:	if set ratios is not a string array, but the path to a file with such an array written line by line
;
; OTIONAL OUTPUTS:
;	index_min:	1-dim index of the best-fit model
;
; COMMON BLOCKS:
;	tt_units:	unit conversion factors
;
; TEST STATUS:
;	checked against by hand calculations
;
; EXAMPLE:
;	fit, '../../tt_data/bptau.ratio','../../tt_data/bptau.dat',/rat,nh=2e21
;
; REQUIREMENTS:
;	In addition to CHINATI PINTofALE is used for the N_H correction.
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 22.09.2005
;	20.09.2007 changed default simdata file to absolute path  HMG
;	20.09.2007 added verbose keyword HMG
;	21.09.2007 ratios is no longer changed on output HMG
;-


pro fit, ratio, datafile,lev=lev,ratiofile=ratiofile,simdata_file=simdata_file, index_min=index_min,c_annotation=c_annotation,nh=nh, verbose=verbose
common tt_units
if n_elements(simdata_file) eq 1 then begin
  restore,simdata_file 
  simdata=data
endif else begin 
  restore,'/data/hspc62/steh305/idl/tt_results/grid/rad6e3newfe.sav'
  simdata=data
endelse  
n_sim=n_elements(simdata.n0)

read_obs_data,datafile,obs,obs_ref,nh=nh
obs.name=strcompress(obs.name,/remove_all)

ratios=ratio
calculate_ratios, ratios,obs,simdata,obsratios,simratios,ratiofile=ratiofile


chi_squared=make_array(n_sim,val=0.)



for i=0, (n_sim-1) do begin

  chi_squared[i]=total(((obsratios[*,0]-simratios[*,i])/obsratios[*,1])^2)
  
endfor

contour, chi_squared,simdata.v0/km,alog10(simdata.n0),/irregular,nlev=5,lev=lev,ytit=textoidl('log_{10}(infall density n_0[cm^{-3}])'), xtit=textoidl('infall velocity v_0 [km/s]'),c_linestyle=[0],c_charsize=.8*!p.charsize,xmargin=[6,2],c_annotation=c_annotation

chi_min=min(chi_squared,index_min,/nan)
print,'Best fit:'
print,format='("v0=",I3," [km/s],  log(n0)=",f6.1," log([1/cm^3])")',simdata.v0[index_min]/km,alog10(simdata.n0[index_min])
print,'with an unreduced chi^2',chi_min

if verbose gt 3 then begin
  print, '    ratio     observed    best-fit model'
  for i=0,n_elements(ratios)-1 do print, format='(a16," & ",f5.2,"\pm",f5.2," & ",f5.2," \\")', ratios[i], obsratios[i,*], simratios[i,index_min]
endif
end