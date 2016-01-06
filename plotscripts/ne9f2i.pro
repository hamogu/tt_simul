pro plot_f2i,element
pmultiold=!p.multi
!p.multi=[0,1,1]
!p.charsize=1.3
restore, '../tt_results/grid/rad8e3.sav'
restore, '../tt_results/grid/rad1e4.sav'
restore, '../tt_results/grid/analysegrid.sav'  ;radtemp=default (6e3)
case element of
  'Ne': begin
	name='Ne IX'
	rad6e3f2i=analysegrid.ne9f2i
	rad8e3f2i=rad8e3.ne9f2i
	rad1e4f2i=rad1e4.ne9f2i
	end
  'O':  begin
  	name='O VII'
	rad6e3f2i=analysegrid.o7f2i
	rad8e3f2i=rad8e3.o7f2i
	rad1e4f2i=rad1e4.o7f2i
	end
  else: begin
  	print, 'Elements not supported'
	exit
	end

contour, rad8e3f2i,alog10(rad8e3.n0),rad8e3.v0/km,/irregular,nlev=5,/follow,lev=[0.1,0.3,0.5,0.75,1.,2.],title=name+' f/i ratios',xtit=textoidl('log_{10}(density [cm^{-3}])'), ytit=textoidl('infall velocity v_0 [km/s]'),c_linestyle=[1]
	contour, rad1e4f2i,alog10(rad1e4.n0),rad1e4.v0/km,/irregular,nlev=5,lev=[0.1,0.3,0.5,0.75,1.,2.],/overplot, c_linestyle=[2]
	contour, rad6e3f2i,alog10(data.n0),data.v0/km,/irregular,nlev=5,lev=[0.1,0.3,0.5,0.75,1.,2.],/overplot, c_linestyle=[0]
	end

legend,['radiation temp','6000 K','8000 K','10000 K'],position=[12.8,300.],linestyle=[-1,0,1,2],charsize=1.2                                           

!p.multi=pmultiold                                            
delvar, pmultiold
end