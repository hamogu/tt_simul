;suited in this way to produce a good .eps outfile
!p.charsize=2.
;plot, min_step,xdata.n_steps,/xlog,ytit='number of steps',xtit='minimum step size [cm]';,xmargin=[20,1]
;plot, min_step,xdata.maxdepth/km,/xlog,ytit='shock depth [km]',xtit='minimum step size [cm]';,xmargin=[20,1]
plot, stepsize,stepdata.n_steps,ytit='number of steps',xtit='maximum relative change';,xmargin=[20,1]
;plot, stepsize,stepdata.olyb/1e3,ytit='intensity of one Fe XVII line !C [arbitrary units]',xtit='maximum relative change', xmargin=[10,3]