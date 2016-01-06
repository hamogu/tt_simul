pro tt_simul_test_start_ioneq;number,output=output
common tt_units
;start with different (and very extreme) ioneqs to check the impotance for the results

;tt_simul,300*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/testneutralexh.ioneq','../tt_results/start_testneutralexh'
;tt_simul,300*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/start_ar'
;tt_simul,300*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/test5.ioneq','../tt_results/start_test5'
tt_simul,300*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/mazzotta45.ioneq','../tt_results/start_mazzotta45'
tt_simul,300*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/mazzotta48.ioneq','../tt_results/start_mazzotta48'
tt_simul,300*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/mazzotta50.ioneq','../tt_results/start_mazzotta50'
end