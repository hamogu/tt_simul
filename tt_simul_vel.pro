pro tt_simul_vel, time=time
common tt_units
common out,out250,out300,out350,out400,out450
time=make_array(7)
time[0]=systime(/seconds)
tt_simul,250*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/vel/250'
time[1]=systime(/seconds)
tt_simul,300*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/vel/300'
time[2]=systime(/seconds)
tt_simul,350*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/vel/350'
time[3]=systime(/seconds)
tt_simul,400*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/vel/400'
time[4]=systime(/seconds)
tt_simul,450*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/vel/450'
time[5]=systime(/seconds)
tt_simul,500*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/vel/500'
time[6]=systime(/seconds)
print, time
end