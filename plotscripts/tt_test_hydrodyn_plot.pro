;suited in this way to produce a good .eps outfile
;set breakpoint in tt_test_hydrodyn and issue these commands
!p.charsize=1.5
plot, x_list/u,analytic,line=0,title='constant electron temperature',ytit='Temperature [K]',xtit='time [sec]',yrange=[0,1.1e7];,xmargin=[20,1]
oplot, x_list/u,t_ion_list,psym=1
oplot, x_list/u,make_array(n_elements(x_list),val=1e7),line=2
legend,['simulation','analytic solution','electron temperature'],position=[.5,3e6],linestyle=[0,0,2],psym=[1,0,0],charsize=1.5
  