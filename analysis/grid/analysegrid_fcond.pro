pro analysegrid_fcond,namelist,data
common tt_constants

num=n_elements(namelist)
data={f:fltarr(num),v0:fltarr(num),n0:fltarr(num)}
;namelist is an array of positions for *.sav files

for i=0,num-1 do begin
  print, 'Working on '+namelist[i]
  restore,namelist[i]
  n=n_elements(hydrodyn[0,*])-1
  ;find middle
  mid=min(abs(hydrodyn[0,*]-hydrodyn[0,0]/2),middle)
  dx=hydrodyn[0,middle]-hydrodyn[0,middle+1]
  dt=hydrodyn[5,middle]-hydrodyn[5,middle+1]
  lambda=9.4+1.5*alog10(hydrodyn[5,middle])-0.5*alog10(hydrodyn[2,middle])
  f_cond=1e-5/lambda*hydrodyn[5,middle]^2.5*dt/dx
  f_flux=0.5*hydrodyn[1,middle]*hydrodyn[3,middle]^3+hydrodyn[2,middle]*hydrodyn[3,middle]*1.5*k_boltz*hydrodyn[5,middle]
  data.f[i]=abs(f_cond)/f_flux
  data.v0[i]=hydrodyn[3,n]
  data.n0[i]=hydrodyn[2,n]
endfor          
end             


; f_cond=1e-5*hydrodyn[4,*]^2.5*dt/dx/(9.4+1.5*alog10(hydrodyn[5,*])-.5*alog10(hydrodyn[2,*]))
;f_kin=hydrodyn[1,*]/2.*hydrodyn[3,*]^3
;f_therm=1.5*k_boltz*hydrodyn[2,*]*hydrodyn[3,*]*hydrodyn[4,*]
