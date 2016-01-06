pro vapec2phot, xcmfile,obsfile

common elements, abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
common tt_constants
common tt_units

d_lambda=0.06

readxcm,xcmfile,xcmfile+'.abund',vapec=vapec,nh=nh

read_obs_data,obsfile,obsdata,obsref,nh=nh

;calculate theoretical coronal contributions for all lines contained in the observed data set
coronalcont=obsdata.flux*0
;test=obsdata.flux*0

;line_intens boils down to mygoft, the alternative useres more standard CHIANTI routines -> easier to maintain in the long run
;ionnames=strarr(n_elements(obsdata.name)-1)
;for i=0,n_elements(obsdata.name)-1 do ionnames[i]=strmid(strtrim(obsdata.name[i],2),0, (strsplit(strtrim(obsdata.name[i],2),' '))[2]-1)
;line_intens, ionnames,obsdata.lambda,vapec[0,*],vapec[1,*]*1e14,d_lambda=0.06

for i = 0, n_elements(obsdata.name)-1 do begin
  spectroscopic2ion,strmid(strtrim(obsdata[i].name,2),0, (strsplit(strtrim(obsdata[i].name,2)+' test',' '))[2]-1), ionname
  convertname,ionname,iz,ion
  ch_synthetic,obsdata[i].lambda-d_lambda/2.,obsdata[i].lambda+d_lambda/2., output=out, sngl_ion=ionname, logt=alog10(vapec[0,*]), logem=alog10(vapec[1,*]*1e14), dens=1e10, ioneq_name='/usr/local/hssoft/chianti/dbase/ioneq/mazzotta_etal_9.ioneq'
  coronalcont[i] = total(out.lines.int)*abund[iz-1]*4.*!pi
  ;ch_synthetic,obsdata[i].lambda-d_lambda/2.,obsdata[i].lambda+d_lambda/2., output=out, sngl_ion=ionname, logt=alog10(2.3e6), logem=alog10(0.00015*1e14), dens=1e10, ioneq_name='/usr/local/hssoft/chianti/dbase/ioneq/mazzotta_etal_9.ioneq'
  ;test[i]=total(out.lines.int)*abund[iz-1]*4.*!pi
endfor

negshockflux=where(coronalcont gt obsdata.flux)
if negshockflux[0] ne -1 then begin
  print,'calculated coronal contribution larger than total flux in'
  print, obsdata[negshockflux].name
endif

print,'  Line       & lambda & obs. flux   & coronal flux \\'

;--- write calculated shock contribution to file + terminal ---
openw,shockfile,obsfile+'shock.dat',/get_lun
for i=0,n_elements(coronalcont)-1 do begin
  phtflux=(obsdata[i].flux-coronalcont[i])/c_light/h_planck*obsdata[i].lambda
  phterr=obsdata[i].error/c_light/h_planck*obsdata[i].lambda
  printf,shockfile, phtflux, phterr, obsdata[i].lambda, obsdata[i].name, format='(2e10.2,f6.2,a'+strtrim(string(strlen(obsdata[i].name)),1)+')'
  print,format='(a12," & ",f6.2," & ",e10.2," & ",e10.2," \\")',strtrim(obsdata[i].name),obsdata[i].lambda,obsdata[i].flux,coronalcont[i]
endfor
printf,shockfile,'-1'
printf,shockfile,'original comment on observed data file:'
printf,shockfile,obsref
printf,shockfile,'produced by Moritz Guenther - vapec2phot'
printf,shockfile,'All observed fluxes have been dereddened with N_H=',nh,'/cm^2'
printf,shockfile,'before substracting the photon flux for a coronal component.'
printf,shockfile,'The coronal component specification is',vapec, ' T in K, norm as in XSPEC'
printf,shockfile,'This file should now contain the shock component only.'
printf,shockfile,'A Note on errors: The vapec model is assumed to be exact, the errors are propagated '
printf,shockfile,'from the line counts only.'
printf,shockfile,'-1'
free_lun,shockfile

end