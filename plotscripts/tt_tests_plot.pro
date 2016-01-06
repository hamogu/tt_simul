!p.charsize=1.8

plot, (findgen(41)+40.)/10.,result[*,16,2],/ylog,yrange=[1e-4,1e0],ytitle='Ionisation fraction',line=2,xmargin=[10,1.5],xrange=[4.,6.], Xtitle='log!D10!N(T in K)', ymargin=[4,1];,title='Ionisation fraction of Cl ions'
oplot, (findgen(41)+40.)/10.,mazzotta[*,16,2],line=0
;oplot, (findgen(41)+40.)/10.,static[*,5,3],line=1
;legend,['Mazzotta et al (1998)','dynamic relaxation','static solution'],psym=[0,0,0],/left,linestyle=[0,2,1],charsize=1.5
for i=3,5 do begin & oplot, (findgen(41)+40.)/10.,result[*,16,i],line=2 & oplot, (findgen(41)+40.)/10.,mazzotta[*,16,i],line=0 & endfor
;for C
; xyouts,5.2,0.005,'C IV'
; xyouts,6.9,0.005,'C VI'
; xyouts,5.3,0.6,'C V'
;legend,['Mazzotta et al','dynamic relax'],psym=[0,0],posit=[6.4,0.6],linestyle=[0,2],charsize=1.5
;legend,['Mazzotta','et al (1998)','dynamic', 'relaxation'],posit=[6.4,0.6],linestyle=[0,-1,2,-1],charsize=1.5
;for Cl
xyouts,4.2,0.2,'Cl III'
xyouts,4.45,0.09,'Cl IV'
xyouts,4.7,0.01,'Cl V'
xyouts,5.05,0.0005,'Cl VI'