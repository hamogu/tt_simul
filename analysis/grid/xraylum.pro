pro xraylum,n0,v0,relxray,maxdepth

common tt_units
namelist=file_search('*.hydrodyn')
for i=0,n_elements(namelist)-1 do namelist[i]=(strsplit(namelist[i],'.',/extract))[0]
n0=make_array(n_elements(namelist))
rho0=n0*0.
v0=n0*0.
maxdepth=n0*0.
relxray=n0*0.

base_abundfile=!xuvtop +'abundance/version_3/grevesse_sauval98.abund'
read_abund,base_abundfile, base_abund,base_abund_ref
read_abund,'tt_simul.abund',abund,abund_ref

min_lambda=1.
max_lambda=50.
ang2kev=12.39854
kev_scale=findgen(12000)*0.00025+0.000125
lambda_scale=ang2kev/kev_scale
lambda_scale=lambda_scale[where(lambda_scale ge min_lambda and lambda_scale le max_lambda)]
n_energies=11008
kev_scale=ang2kev/lambda_scale
elements=['C','N','O','Ne','Mg','Si','S','Fe']
elemnumb=[6,7,8,10,12,14,16,26]-1

for ii=0,n_elements(namelist)-1 do begin
  restore, namelist[ii]+'.sav'
  rho0[ii]=hydrodyn[1,n_elements(hydrodyn[0,*])-1]
  n0[ii]=alog10(hydrodyn[2,n_elements(hydrodyn[0,*])-1])
  v0[ii]=hydrodyn[3,n_elements(hydrodyn[0,*])-1]/km
  maxdepth[ii]=hydrodyn[0,0]/km
  restore, namelist[ii]+'spectrum.sav'
  index=where(base_spectra.intpspec eq base_spectra.intpspec)
  xrayloss=total(kev_scale[index]*base_spectra.intpspec[index])
  print, n0[ii],v0[ii],rho0[ii]
  print, xrayloss
  for j=0,7 do begin
    index=where(base_spectra.(j+2) eq base_spectra.(j+2))
    loss=base_abund[elemnumb[j]]/abund[elemnumb[j]]*total(kev_scale*(base_spectra.(j+2))[index])
    print, loss
    xrayloss+=loss
  endfor
  print, xrayloss*1.602e-9, 0.5*(v0[ii]*km)^3.*rho0[ii]*0.75
  relxray[ii]=xrayloss*1.602e-9/(0.5*(v0[ii]*km)^3.*rho0[ii]*0.75)
endfor

end