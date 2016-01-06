pro analysegrid_depth,namelist,data
common tt_constants

num=n_elements(namelist)
data={maxdepth:fltarr(num),v0:fltarr(num),n0:fltarr(num)}
;namelist is an array of positions for *.sav files

for i=0,num-1 do begin
  print, 'Working on '+namelist[i]
  if file_test('./deep/'+namelist[i]) then restore,'./deep/'+namelist[i] else restore, namelist[i]
  n=n_elements(hydrodyn[0,*])-1
  ;find middle
  data.maxdepth[i]=hydrodyn[0,0]
  data.v0[i]=hydrodyn[3,n]/1e5
  data.n0[i]=alog10(hydrodyn[2,n])
endfor          
end             


; f_cond=1e-5*hydrodyn[4,*]^2.5*dt/dx/(9.4+1.5*alog10(hydrodyn[5,*])-.5*alog10(hydrodyn[2,*]))
;f_kin=hydrodyn[1,*]/2.*hydrodyn[3,*]^3
;f_therm=1.5*k_boltz*hydrodyn[2,*]*hydrodyn[3,*]*hydrodyn[4,*]
