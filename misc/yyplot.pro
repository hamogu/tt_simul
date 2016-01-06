;+
; NAME:
;	yyplot
;
; PURPOSE:
;	x-y-plot for 2 data sets with different y axises
;
; CATEGORY:
;	plot tool
;
; CALLING SEQUENCE:
;	yyplot, x,yl,yr,xlog=xlog,ylog=ylog,xtitle=xtitle,yltitle=yltitle,yrtitle=yrtitle,xrange=xrange,ylreange=ylrange,yrrange=yrrange
;
; INPUTS:
;	x:	data for xaxis
;	yl:	data with left y-axis
;	yr:	same for right
;
; OPTIONAL INPUTS:
;	some of the usual graphics keywords with small modification
;	yrange:	here yrrange for right axis, ylrange for left
;	ytitle: analog, yltitle and yrtitle
;	
; KEYWORD PARAMETERS:
;	/xlog, /ylog
;
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther	17-05-2005
;-

pro yyplot, x,yl,yr,xlog=xlog,ylog=ylog,xtitle=xtitle,yltitle=yltitle,yrtitle=yrtitle,xrange=xrange,ylrange=ylrange,yrrange=yrrange,$
	psyml=psyml,psymr=psymr,yllog=yllog,yrlog=yrlog,xmargin=xmargin,ymargin=ymargin,xstyle=xstyle,ystyle=ystyle,noerase=noerase,$
	position=position ,charsize=charsize,lcharsize=lcharsize,rcharsize=rcharsize,xtickv=xtickv,xcharsize=xcharsize,CHARTHICK=CHARTHICK,$
	THICK=THICK,XCOLOR=XCOLOR,yCOLOR=yCOLOR,linethick=linethick
  if keyword_set(ylog) then begin
    yrlog=ylog
    yllog=ylog
  endif 
  if n_elements(xmargin) eq 0 then xmargin=[8,8]
  if n_elements(ymargin) eq 0 then ymargin=[4,4]
  if n_elements(ystyle) eq 0 then ystyle=0
  if n_elements(lcharsize) eq 0 and n_elements(charsize) ne 0 then lcharsize=charsize
  if n_elements(rcharsize) eq 0 and n_elements(charsize) ne 0 then rcharsize=charsize
  plot, x,yl,xlog=xlog,ylog=yllog,xtitle=xtitle,ytitle=yltitle,xrange=xrange,yrange=ylrange, ystyle=8+ystyle,psym=psyml,xmargin=xmargin,ymargin=ymargin,xstyle=xstyle, noerase=noerase, position=position, charsize=lcharsize, xtickv=xtickv, xcharsize=xcharsize, xthick=linethick,ythick=linethick,thick=linethick,charthick=linethick
  if n_elements(yrrange) eq 0 then yrrange=[min(yr),max(yr)]
  if n_elements(ycolor) ne 0 then begin
    yline=ycolor eq xcolor ? 2 : 0 
  endif
  axis,yaxis=1,yrange=yrrange,ytitle=yrtitle,ystyle=1+ystyle,/save,ylog=yrlog,xmargin=xmargin,ymargin=ymargin,charsize=rcharsize,xtickv=xtickv,xcharsize=xcharsize,color=ycolor,ythick=linethick,charthick=linethick
  
  oplot, x,yr,line=yline,psym=psymr,color=ycolor,thick=linethick
end

