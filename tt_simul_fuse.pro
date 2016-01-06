pro tt_simul_fuse, time=time
common tt_units
time=make_array(7)
time[0]=systime(/seconds)
;tt_simul,525*km/s,1.67e-11*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/twhya_uv/525',abund='../tt_data/twhya2.abund'
;time[1]=systime(/seconds)
;tt_simul,450*km/s,1.1e-11*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/twhya_uv/450',abund='../tt_data/twhya2.abund'
;time[2]=systime(/seconds)
tt_simul,350*km/s,5.5e-12*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/twhya_uv/350',abund='../tt_data/twhya2.abund'
time[3]=systime(/seconds)
tt_simul,250*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/twhya_uv/250',abund='../tt_data/twhya2.abund'
time[4]=systime(/seconds)
print, time
end
