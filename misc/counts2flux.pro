pro counts2flux,datafile,corafile


read_obs_data,datafile+'.counts',obs_data,ref,/no_processing

cora=read_ascii(corafile,data_start=1,header=dummy)
exptime=float(strmid(dummy,2))
print,'Exposure time [in s]:',exptime
cora=cora.field1

;write new dat file with line fluxes instead of line counts
openw,lun,datafile+'.dat',/get_lun
for i=0,n_elements(obs_data.flux)-1 do begin
  dummy=min(abs(cora[0,*]-obs_data[i].lambda),index)
  printf,lun,format='(2e8.1,f6.2,a)', obs_data[i].flux/exptime/cora[5,index], obs_data[i].error/exptime/cora[5,index], obs_data[i].lambda, obs_data[i].name
endfor

index=where(1 eq strmatch(ref,'*unit*'))
ref[index]='units: cts/s/cm^2'
printf,lun,'-1'
for i=0,n_elements(ref)-1 do begin
  printf,lun,ref[i]
endfor
printf,lun,'-1'
free_lun,lun

end