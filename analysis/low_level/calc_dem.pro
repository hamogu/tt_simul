function calc_dem, filename, nbins,log=log

restore, filename

prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens
n_steps=n_elements(hydrodyn[5,*])
volem=volem[2:n_steps-2] ;last step includes temP2e4, where simulation is not valid
if ~keyword_set(log) then hist1=histogram(hydrodyn[5,2:n_steps-2],nbins=nbins,reverse=R,locations=Temp) else hist1=histogram(alog10(hydrodyn[5,1:n_steps-2]),nbins=nbins,reverse=R,locations=Temp)

dem=make_array(nbins)
for i=0,nbins-1 do dem[i]= r[i] ne r[i+1] ? total(volem[R[R[i]:R[i+1]-1]]):0
return,[[temp],[dem]]
end