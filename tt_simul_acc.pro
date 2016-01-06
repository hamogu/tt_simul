pro tt_simul_acc, paramfile,time=time
common tt_units
time=make_array(n_elements(paramfile)+1)
time[0]=systime(/seconds)
for i=0,n_elements(paramfile)-1 do begin
  tt_simul,300*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/acc/'+file_basename(paramfile[i],'.par'),paramfile=paramfile[i]
  time[i+1]=systime(/seconds)
endfor
print, time
end