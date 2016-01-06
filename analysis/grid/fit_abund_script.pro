index_min=61 ;from fit.pro
fit_abund, '../../tt_data/ne2o.ratio','../../tt_data/twhya.dat',index_min,/ratiofile
fit_abund, '../../tt_data/n2o.ratio','../../tt_data/twhya.dat',index_min,/ratiofile
fit_abund, ['clya/olya','clya/o7r','clya/o7i','clya/o7f'],'../../tt_data/twhya.dat',index_min
fit_abund, ['fe17l1705/olya','fe17l1705/o7r','fe17l1705/o7i','fe17l1705/o7f'],'../../tt_data/twhya.dat',index_min
