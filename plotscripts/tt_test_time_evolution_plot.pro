!p.charsize=1.7
plot, findgen(41)/10.,make_array(41,val=1.-n_inf),line=2,xtitle='Time [timescales]',ytitle='n!DHI!N/n!DH!N' ;,title='time development of ionisation stages of H'
oplot,findgen(41)/10.,nh1,psym=1
t=findgen(41)/10.
oplot,t , 1.-n_inf/(1.-(1.-n_inf)*exp(-t)),line=0
legend,['simulation','analytic solution','equilibrium value'],psym=[1,0,0],position=[0.8,0.2],linestyle=[0,0,2],charsize=1.7
