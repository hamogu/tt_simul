pro plot_lineratio,ratio,title=title,lev=lev,old=old,legendpos=legendpos,c_annotation=c_annotation,xmargin=xmargin
common tt_units
pmultiold=!p.multi
!p.multi=[0,1,1]
; !p.charsize=1.3
if keyword_set(old) then begin
  restore, '../tt_results/grid/rad8e3.sav'
  restore, '../tt_results/grid/rad1e4.sav'
  restore, '../tt_results/grid/rad6e3.sav'
endif else begin
  restore, '../tt_results/grid/rad8e3newfe.sav'
  rad8e3=data
  restore, '../tt_results/grid/rad10e3newfe.sav'
  rad1e4=data
  restore, '../tt_results/grid/rad6e3newfe.sav'
  rad6e3=data
endelse   
;restore, '../tt_results/grid/analysegrid.sav'  ;radtemp=default (6e3)
radratio=make_array(3,n_elements(rad6e3.o7f2i))
datastructure=['rad6e3','rad8e3','rad1e4']

case strlowcase(ratio) of
  'ne9f2i': begin
	name='Ne IX f/i'
	radratio[0,*]=rad6e3.ne9f2i
	radratio[1,*]=rad8e3.ne9f2i
	radratio[2,*]=rad1e4.ne9f2i
	lev=[0.1,0.3,0.5,0.75,1.,2.]
	end
  'o7f2i':  begin
  	name='O VII f/i'
	radratio[0,*]=rad6e3.o7f2i
	radratio[1,*]=rad8e3.o7f2i
	radratio[2,*]=rad1e4.o7f2i
	lev=[0.05,0.1,0.3,0.5,0.75,1.]
	end
  'o7g':  begin
  	name='O VII (f+i)/r'
	radratio[0,*]=(rad6e3.o7f+rad6e3.o7i)/rad6e3.o7r
	radratio[1,*]=(rad8e3.o7f+rad8e3.o7i)/rad8e3.o7r
	radratio[2,*]=(rad1e4.o7f+rad1e4.o7i)/rad1e4.o7r
	follow=follow
	lev=[0.6,0.8,1.,1.2]
	end
  'ne9g':  begin
  	name='Ne IX (f+i/r)'
	radratio[0,*]=(rad6e3.ne9f+rad6e3.ne9i)/rad6e3.ne9r
	radratio[1,*]=(rad8e3.ne9f+rad8e3.ne9i)/rad8e3.ne9r
	radratio[2,*]=(rad1e4.ne9f+rad1e4.ne9i)/rad1e4.ne9r
	end
  'ne9gblend':  begin
  	name='Ne IX (f+i/r) incl. the blended Fe XIX 13.52'
	radratio[0,*]=(rad6e3.ne9f+rad6e3.ne9i+rad6e3.fe19l1352)/rad6e3.ne9r
	radratio[1,*]=(rad8e3.ne9f+rad8e3.ne9i+rad8e3.fe19l1352)/rad8e3.ne9r
	radratio[2,*]=(rad1e4.ne9f+rad1e4.ne9i+rad1e4.fe19l1352)/rad1e4.ne9r
	end
	
  'ne9f2iblend' : begin
  	name='Ne IX f/i incl. the blended Fe XIX 13.52'
	radratio[0,*]=rad6e3.ne9f/(rad6e3.ne9i+rad6e3.fe19l1352)
	radratio[1,*]=rad8e3.ne9f/(rad8e3.ne9i+rad8e3.fe19l1352)
	radratio[2,*]=rad1e4.ne9f/(rad1e4.ne9i+rad1e4.fe19l1352)
	end
  
;  'nelya2olya': begin
;	name=textoidl('Ne Lyman \alpha / O Lyman \alpha')
;	radratio[0]=data.nelya/data.olya
;	radratio[1]=rad8e3.nelya/rad8e3.olya
;	radratio[2]=rad1e4.nelya/rad1e4.olya
;	end
	
  else: begin
  	;make this procedure flexible, but it is a waste of time to sort all possible combinations into a proper title, use title keyword
	names=strsplit(ratio,'/',/extract)
	test=execute('tagnames=strlowcase(tag_names('+datastructure[0]+'))' )
	if ((where(tagnames eq strlowcase(names[0])) ne -1) and (where(tagnames eq strlowcase(names[1])) ne -1)) then begin
	  for i=0,n_elements(radratio[*,0])-1 do test=execute('radratio[i,*]='+datastructure[i]+'.'+names[0]+'/' +datastructure[i]+'.'+names[1])
	  if ratio eq 'nelya/olya' then title=textoidl('Ne Lyman \alpha / O Lyman \alpha')
	  name= n_elements(title) ne 0 ? title : ratio & print, 'Use optional input title= to set the title'
	endif else begin
	  print, 'Ratio not supported - Tag names do not exist'
	  return
	endelse  
	end
endcase

contour, radratio[0,*],rad8e3.v0/km,alog10(rad6e3.n0),/irregular,nlev=5,/follow,lev=lev,title=name+' ratios',ytit=textoidl('log_{10}(infall density n_0[cm^{-3}])'), xtit=textoidl('infall velocity v_0 [km/s]'),c_linestyle=[0],c_charsize=1.5,c_annotation=c_annotation,xmargin=xmargin

;n0=alog10(rad6e3.n0)
;v0=rad6e3.v0/km
;contour, sort2d(radratio[1,*],v0,n0),sort2d(v0,v0,n0),sort2d(n0,v0,n0),lev=[.95,1.25],/noerase,/cell_fill, c_line=[1,-1],c_orient=[45,0],c_thick=[1,1e-8],c_spac=[.15,100]


contour, radratio[2,*],rad1e4.v0/km,alog10(rad1e4.n0),/irregular,nlev=5,lev=lev,/overplot, c_linestyle=[2];,follow=follow,c_charsize=1.3

contour, radratio[1,*],rad6e3.v0/km,alog10(rad8e3.n0),/irregular,nlev=5,lev=lev,/overplot, c_linestyle=[1];,follow=follow,c_charsize=1.3	



;legend,['radiation temp','6000 K','8000 K','10000 K'],position=legendpos,linestyle=[-1,0,1,2],charsize=!p.charsize*0.85
print, 'Plotted are ENERGY flux ratios'
!p.multi=pmultiold                                            

end