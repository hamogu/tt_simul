pro analyse_hydrodyn,hydrodyn
common tt_units
!p.multi=[0,2,2]
plot, hydrodyn[0,*],hydrodyn[4,*],title='ion and electron temperature'
oplot, hydrodyn[0,*],hydrodyn[5,*],line=2
plot, hydrodyn[0,*],hydrodyn[4,*],title='ion and electron temperature',/xlog,/ylog,xrange=[1e1,1e8]
oplot, hydrodyn[0,*],hydrodyn[5,*],line=2
plot, hydrodyn[0,*],hydrodyn[2,*],title='density',/ylog
plot, hydrodyn[0,*],hydrodyn[4,*],title='ion and electron temperature',xrange=[-1e0,1e5]
oplot, hydrodyn[0,*],hydrodyn[5,*],line=2
print,'End depth:', max(hydrodyn[0,*])/km,' km'
print, 'Maximum temperature in K: ions	   ',max(hydrodyn[4,*])
print, '			 electrons:',max(hydrodyn[5,*])
print,'number of steps',n_elements(hydrodyn[0,*])
index=where(hydrodyn[5,*] ne hydrodyn[5,*],ct)
print,'errors rows in hydrodyn (NaN)',ct
end