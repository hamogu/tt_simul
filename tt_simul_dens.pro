pro tt_simul_dens, time=time
common tt_units
time=make_array(7)
time[0]=systime(/seconds)
tt_simul,300*km/s,1.e-10*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/dens/1e10'
time[1]=systime(/seconds)
tt_simul,300*km/s,1.e-11*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/dens/1e11'
time[2]=systime(/seconds)
tt_simul,300*km/s,1.e-12*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/dens/1e12'
time[3]=systime(/seconds)
tt_simul,300*km/s,1.e-13*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/dens/1e13'
time[4]=systime(/seconds)
tt_simul,300*km/s,1.e-14*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/dens/1e14'
time[5]=systime(/seconds)
tt_simul,300*km/s,5.e-13*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/dens/5e13'
time[6]=systime(/seconds)
print, time
end