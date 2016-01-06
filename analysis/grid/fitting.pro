pro fitting, lev=lev, bptau=bptau,kastner=kastner, old=old
common tt_units
if keyword_set(old) then begin
  restore,'rad6e3.sav' 
endif else begin 
  restore,'rad6e3newfe.sav'
  rad6e3=data
endelse  
;try fitting of TW Hya line ratios to my model
;simple by hand start
;a[0] is line, a[1] is error
;---TW HYA
olya=[30.4,2.8]/18.97	; flux/wavelength = energy
o7r=[17.4,2.1]/21.61
o7i=[9.5,1.5]/21.80
o7f=[0.4,0.6]/22.1
nelya=[6.3,1.5]/12.13
ne9r=[23.1,2.2]/13.45
ne9i=[17.4,1.1]/13.55
ne9f=[8.2,1.5]/13.70
if keyword_set(bptau) then begin
;---BP TAU
olya=[14.2,1.7]/18.97	; flux/wavelength = energy
o7r=[8.9,1.6]/21.61
o7i=[7.1,1.5]/21.80
o7f=[2.7,1.0]/22.1
nelya=[13.2,1.7]/12.13
ne9r=[5.2,1.4]/13.45
ne9i=[4.4,1.4]/13.55
ne9f=[1.7,1.0]/13.70
endif
if keyword_set(kastner) then begin
;---TW Hya nach Kastner et al (2002)
olya=[195.7,28.0]/18.97	;count/wavelength = energy
o7r=[91.3,36.3]/21.61
o7i=[105.4,40.4]/21.80
o7f=[01.0,20.5]/22.1   ;estimated, article just gives a 2 sigma upper limit, but I cannot run with an 0 input
nelya=[73.7,7.8]/12.13
ne9r=[122.4,11.9]/13.45
ne9i=[79.2,10.5]/13.55
ne9f=[35.3,8.6]/13.70
endif
o7f2i=[o7f[0]/o7i[0],sqrt(o7f[1]^2/o7i[0]^2+o7f[0]^2/o7i[0]^4*o7i[1]^2)]
o7g=[(o7f[0]+o7i[0])/o7r[0],sqrt((o7f[1]^2+o7i[1]^2)/o7r[0]^2+(o7f[0]+o7i[0])^2/o7r[0]^4*o7r[1]^2)]
olya2o7r=[olya[0]/o7r[0],sqrt(olya[1]^2/o7r[0]^2+olya[0]^2/o7r[0]^4*o7r[1]^2)]


ne9f2i=[ne9f[0]/ne9i[0],sqrt(ne9f[1]^2/ne9i[0]^2+ne9f[0]^2/ne9i[0]^4*ne9i[1]^2)]
ne9g=[(ne9f[0]+ne9i[0])/ne9r[0],sqrt((ne9f[1]^2+ne9i[1]^2)/ne9r[0]^2+(ne9f[0]+ne9i[0])^2/ne9r[0]^4*ne9r[1]^2)]
nelya2ne9r=[nelya[0]/ne9r[0],sqrt(nelya[1]^2/ne9r[0]^2+nelya[0]^2/ne9r[0]^4*ne9r[1]^2)]

xi_squared=make_array(n_elements(rad6e3.o7i),val=0.)

for i=0, (n_elements(rad6e3.o7i)-1) do begin
  xi_squared[i]+=(o7f2i[0]-rad6e3.o7f2i[i])^2/o7f2i[1]^2
  xi_squared[i]+=(o7g[0]-(rad6e3.o7f[i]+rad6e3.o7i[i])/rad6e3.o7r[i])^2/o7g[1]^2
  xi_squared[i]+=(olya2o7r[0]-(rad6e3.olya[i]/rad6e3.o7r[i]))^2/olya2o7r[1]^2
  
  xi_squared[i]+=(ne9f2i[0]-rad6e3.ne9f2i[i])^2/ne9f2i[1]^2
  xi_squared[i]+=(ne9g[0]-(rad6e3.ne9f[i]+rad6e3.ne9i[i])/rad6e3.ne9r[i])^2/ne9g[1]^2
  xi_squared[i]+=(nelya2ne9r[0]-(rad6e3.nelya[i]/rad6e3.ne9r[i]))^2/nelya2ne9r[1]^2
endfor

contour, xi_squared,alog10(rad6e3.n0),rad6e3.v0/km,/irregular,nlev=5,/follow,lev=lev,xtit=textoidl('log_{10}(infall density n_0[cm^{-3}])'), ytit=textoidl('infall velocity v_0 [km/s]'),c_linestyle=[1],c_charsize=1.0

xi_min=min(xi_squared,index_min,/nan)
print,'Best fit:'
print,'v0=',rad6e3.v0[index_min],'    n0=',rad6e3.n0[index_min]
print,'with an unreduced xi^2',xi_min

end