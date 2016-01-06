;suited in this way to produce a good .eps outfile
!p.charsize=1.8
plot, hydrodyn005[0,*]/km,ioneq_005[*,25,16],xtit='depth [km]',ytit='ionisation fraction of Fe XVII'
oplot, hydrodyn01[0,*]/km,ioneq_01[*,25,16],line=1
oplot, hydrodyn04[0,*]/km,ioneq_04[*,25,16],line=2
oplot, hydrodyn07[0,*]/km,ioneq_07[*,25,16],line=3
legend,['max relative change per step','5 %','10 %','40 %','70 %'],/right,linestyle=[-1,0,1,2,3],charsize=1.4,box=0                                            