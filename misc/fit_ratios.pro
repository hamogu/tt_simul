;The quadratic form H ist distributed as chi^2 (for the usual assumption of Gaussian errors etc.)
; So in principle this method not only deliveres the best fit, but also error ranges and the measure of the goodness-of-fit
;Unfortunately the covariance matrix tends to be numerically difficult to invert and despite trying invert and la_invert
; with the base matrix v and with 10.*v or 100.*v to make the numbers closer to one
;I could not get a sufficiently accurate invers v, numerical errors always result in h gt 1e9!!!
;-> abolish this method and to a simpler quadratic form!
function quad_form_h, flux,error,mod_flux

;--- ratind contains the ratio_p=flux_i/flux_j, where, p,i,j are stored in the matirx ratind[p-1,*]=[p,i,j]
if n_elements(flux) ne n_elements(error) then message,'Dimensions of flux and error do not agree!'
n=n_elements(flux)
if n gt 1 then begin
  ratind=make_array(n*(n-1)/2,3,/int)
  ratind[*,0]=indgen(n*(n-1)/2)
  p=0
  for i= 1, n-1 do for j= i+1, n do begin & ratind[p,*]=[p+1,i,j] & p++ & end

  obs_ratio=alog(flux[ratind[*,1]-1]/flux[ratind[*,2]-1])
  mod_ratio=alog(mod_flux[ratind[*,1]-1]/mod_flux[ratind[*,2]-1])
  sig2=error^2./flux^2.

  ;--- calculate covariance matrix
  V=make_array(n*(n-1)/2,n*(n-1)/2)
  for i=0, n*(n-1)/2-1 do for j=0, n*(n-1)/2-1 do v[i,j]= (ratind[i,1] eq ratind[j,1]) * sig2[ratind[i,1]-1] +  (ratind[i,2] eq ratind[j,2]) * sig2[ratind[i,2]-1] -  (ratind[i,1] eq ratind[j,2]) * sig2[ratind[i,1]-1] -  (ratind[i,2] eq ratind[j,1]) * sig2[ratind[i,2]-1]

  ;H=(obs_ratio-model_ratio) ## invert(v,/double) ## transpose(obs_ratio-model_ratio)
  ;make use of the matrix_multiply which is recommended instead of ##  in the IDL manual
  H= matrix_multiply(matrix_multiply((obs_ratio-mod_ratio),invert(v,/double),/atranspose),(obs_ratio-mod_ratio))
  return, H
 endif else return,0.
end


function quad_form, flux,error,mod_flux

;--- ratind contains the ratio_p=flux_i/flux_j, where, p,i,j are stored in the matirx ratind[p-1,*]=[p,i,j]
if n_elements(flux) ne n_elements(error) then message,'Dimensions of flux and error do not agree!'
n=n_elements(flux)
if n gt 1 then begin
  ratind=make_array(n*(n-1)/2,3,/int)
  ratind[*,0]=indgen(n*(n-1)/2)
  p=0
  for i= 1, n-1 do for j= i+1, n do begin & ratind[p,*]=[p+1,i,j] & p++ & end

  obs_ratio=alog(flux[ratind[*,1]-1]/flux[ratind[*,2]-1])
  mod_ratio=alog(mod_flux[ratind[*,1]-1]/mod_flux[ratind[*,2]-1])
  sig2=error^2./flux^2.
  return, total((obs_ratio-mod_ratio)^2/sig2)
 endif else return,0.
end


; unfertig, fehlt noch die Beschaffung der model grids -> muss erst mal das programmieren!!
pro fit_ratios, obsfile,model,lev=lev,c_annotation=c_annotation
common tt_units

n_models=n_elements(model->n0())
h=make_array(n_models)

read_obs_data, obsfile,obsdata,ref
index=where(1 eq strmatch(ref,'*Observation:*'))
if index[0] ne -1 then begin
  print,ref[index]
  if strmatch(ref[index],'*XMM*') then dlambda=0.06
  if strmatch(ref[index],'*Chandra/HEG*') then dlambda=0.01
  if strmatch(ref[index],'*Chandra/MEG*') then dlambda=0.02
  print,'Setting dlambda to ',dlambda
endif else begin
  dlambda=0.06
  print, 'No information about instrument found in obsfile, setting dlambda=0.06 (XMM RGS)'
endelse
obsdata=obsdata[SORT(obsdata.iz)]
izelem=[UNIQ(obsdata.iz)]
n_elem=n_elements(izelem)
for i=0,n_elem-1 do begin  ;loop over elements 
  print,'Working on element ', obsdata[izelem[i]].iz
  index=where(obsdata.iz eq obsdata[izelem[i]].iz)
  model_flux=make_array(n_elements(index))
  for j=0,n_models-1 do begin  ;loop over models
    if n_elements(index) gt 1 then begin
      for k=0,n_elements(index)-1 do model_flux[k]= model->getlineflux(obsdata[index[k]].iz, obsdata[index[k]].ion, obsdata[index[k]].lambda, dlambda,j) ;loop over lines
      h[j]+=quad_form(obsdata[index].flux,obsdata[index].error,model_flux)
    end ;else there is no information (one line only), so we can save the time to calculate that line
  endfor
endfor

h_red=h/(n_elements(obsdata.lambda)-2)
contour, h_red,model->v0()/km,alog10(model->n0()),/irregular,nlev=5,lev=lev,ytit=textoidl('log_{10}(infall density n_0[cm^{-3}])'), xtit=textoidl('infall velocity v_0 [km/s]'),c_linestyle=[0],c_charsize=.8*!p.charsize,xmargin=[6,2],c_annotation=c_annotation

chi_min=min(h_red,index_min,/nan)
print,'Best fit:'
print,format='("v0=",I3," [km/s],  log(n0)=",f6.1," log([1/cm^3])")', (model->v0())[index_min]/km,alog10((model->n0())[index_min])
print,'with lowest quad. form',chi_min

end