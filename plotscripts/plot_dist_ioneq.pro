!p.charsize=1.8
;step 0 is the last step, not the first
stepnumber=376-findgen(376)
;yyplot, stepnumber,dist_ioneq(ioneqar[0:375,*,*],reform(hydar[5,0:375]))/30./2,delta_x[0:375]/km,/yrlog,yrrange=[8e-4,5.],Yltitle='difference to equilibrium (line)',yrtitle='step size [km] (dashed)',xtit='number of steps',xmargin=[9,9],xrange=[0,375],xstyle=1,ymargin=[4,1]
yyplot, hydar[0,0:375]/km,dist_ioneq(ioneqar[0:375,*,*],reform(hydar[5,0:375]))/30./2,delta_x[0:375]/km,/yrlog,yrrange=[8e-4,5.],Yltitle='difference to equilibrium (line)',yrtitle='step size [km] (dashed)',xtit='depth [km]',xmargin=[9,9],ymargin=[4,1]
