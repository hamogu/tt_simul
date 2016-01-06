pro tt_simul_run,number,output=output
common tt_units
tt_simul,300*km/s,1.67e-12*g/cm^3,2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/run'+strtrim(string(number),2),output=output
end

;-1.00000e+00 2.08316e-12 1.00759e+12 5.25000e+07 2.00000e+04 2.00000e+04 8.60519e-01 0.00000e+00
.comp /data/hspc62/steh305/idl/tt_simul/analysis/low_level/correct_pops.pro
!XUVTOP = '/usr/local/hssoft/chianti_5.2/dbase/'

tt_simul,525*km/s,2.08316e-12*g/cm^3,2e4,2e4,'../../tt_data/arnaud_rothenflug_ext43.ioneq','.',output=output, abund = 'tt_simul.abund',/quiet