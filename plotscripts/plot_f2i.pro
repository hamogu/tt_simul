pro plot_lineratio,ratio
common tt_units
pmultiold=!p.multi
!p.multi=[0,1,1]
!p.charsize=1.3
restore, '../tt_results/grid/rad8e3.sav'
restore, '../tt_results/grid/rad1e4.sav'
restore, '../tt_results/grid/analysegrid.sav'  ;radtemp=default (6e3)
case strlowcase(element) of
  'ne9f2i': begin
	name='Ne IX f/i'
	rad6e3ratio=data.ne9f2i
	rad8e3ratio=rad8e3.ne9f2i
	rad1e4ratio=rad1e4.ne9f2i
	end
  'o7f2i':  begin
  	name='O VII f/i'
	rad6e3ratio=data.o7f2i
	rad8e3ratio=rad8e3.o7f2i
	rad1e4ratio=rad1e4.o7f2i
	end
  'o7g':  begin
  	name='O VII (f+i)/r'
	rad6e3ratio=(data.o7f+data.o7i)/data.o7r
	rad8e3ratio=(rad8e3.o7f+rad8e3.o7i)/rad8e3.o7r
  	rad1e4ratio=(rad1e4.o7f+rad1e4.o7i)/rad1e4.o7r
	end
  'ne9g':  begin
  	name='Ne IX (f+i/r)'
	rad6e3ratio=(data.ne9f+data.ne9i)/data.ne9r
	rad8e3ratio=(rad8e3.ne9f+rad8e3.ne9i)/rad8e3.ne9r
  	rad1e4ratio=(rad1e4.ne9f+rad1e4.ne9i)/rad1e4.ne9r
	end
  else: begin
  	print, 'Ratio not supported'
	return
	end
endcase

contour, rad8e3ratio,alog10(rad8e3.n0),rad8e3.v0/km,/irregular,nlev=5,/follow,lev=[0.1,0.3,0.5,0.75,1.,2.],title=name+' f/i ratios',xtit=textoidl('log_{10}(infall density n_0[cm^{-3}])'), ytit=textoidl('infall velocity v_0 [km/s]'),c_linestyle=[1]
	contour, rad1e4ratio,alog10(rad1e4.n0),rad1e4.v0/km,/irregular,nlev=5,lev=[0.1,0.3,0.5,0.75,1.,2.],/overplot, c_linestyle=[2],/follow
	contour, rad6e3ration,alog10(data.n0),data.v0/km,/irregular,nlev=5,lev=[0.1,0.3,0.5,0.75,1.,2.],/overplot, c_linestyle=[0],/follow
	

legend,['radiation temp','6000 K','8000 K','10000 K'],position=[12.8,300.],linestyle=[-1,0,1,2],charsize=1.2                                           

!p.multi=pmultiold                                            

end