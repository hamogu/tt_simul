pro grid_point::newline, iz,ion,lambda,flux
  if n_params() eq 3 then begin
    common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
    restore,self.file
    read_ioneq, strmid(self.file,0,strlen(self.file)-4)+'.ioneq' ,ioneq_logt,ioneq,ioneq_ref
    prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens
    n=n_elements(hydrodyn[0,*])-1
    flux=total(volem*calcmygoft(iz,ion,[lambda-.01,lambda+.01],reform(elecdens),reform(hydrodyn[5,*]),/quiet))
    zion2spectroscopic,iz,ion,ionname
    print,'Calculating ',ionname,' line at ',lambda, 'for n0=',self.n0,' and v0=',self.v0 
  endif
  newline={line,iz:iz,ion:ion,lambda:lambda,flux:flux}
  ;print, 'Adding line to n0=',self.n0,' and v0=',self.v0, ': ', newline
  self.plines=ptr_valid(self.plines) ? ptr_new([*self.plines,newline]) : ptr_new(newline)
  self.saved=0
end

function grid_point::getlineflux, iz,ion,lambda,dlambda
  if ptr_valid(self.plines) then index=where((*self.plines).iz eq iz and (*self.plines).ion eq ion and (*self.plines).lambda le lambda+dlambda and (*self.plines).lambda ge lambda-dlambda) else index=[-1]
  if index[0] eq -1 then begin
    self->newline,iz,ion,lambda
    index=n_elements(*self.plines)-1
  end
  return, total((*self.plines)[index].flux)
end

function grid_point::saved
  return, self.saved
end

function grid_point::n0
  return, self.n0
end

function grid_point::v0
  return, self.v0
end

pro grid_point::sethydparam,Te_0,Tion_0,maxdepth,Te_max,Tion_max,n_steps, diff2eq
  self.Te_0=Te_0
  self.Tion_0=Tion_0
  self.maxdepth=maxdepth
  self.Te_max=Te_max
  self.Tion_max=Tion_max
  self.n_steps=n_steps
  self.diff2eq=diff2eq
  self.saved=0
end

pro grid_point::gethydparam,Te_0,Tion_0,maxdepth,Te_max,Tion_max,n_steps, diff2eq
  Te_0=self.Te_0
  Tion_0=self.Tion_0
  maxdepth=self.maxdepth
  Te_max=self.Te_max
  Tion_max=self.Tion_max
  n_steps=self.n_steps
  diff2eq=self.diff2eq
end

pro grid_point::CLEANUP
  IF (PTR_VALID(self.plines)) THEN PTR_FREE,self.plines
end

function grid_point::INIT,n0,v0,file
  self.n0=n0
  self.v0=v0
  if n_params() eq 2 then begin
   self.file=v0n02gridfile(v0,n0)+'.sav'
   message,'file: '+self.file+' for (v0,n0):'+string(v0)+' , '+string(n0),/info
  endif else self.file=file
  self.saved=0
  return,1
end

pro grid_point__define
  A={line,iz:0,ion:0,lambda:0.,flux:0.}
  void={grid_point,v0:0.,n0:0.,Te_0:0.,Tion_0:0.,maxdepth:0.,Te_max:0.,Tion_max:0.,n_steps:0,diff2eq:0.,file:'', saved:1,plines:ptr_new()}
end

function model_grid::n0
  if ptr_valid(self.pgp) then begin
    n=n_elements(*self.pgp)
    result=make_array(n)
    for i=0,n-1 do result[i]=(*self.pgp)[i]->n0()
    return, result
  end else return, -1
end

function model_grid::v0
  if ptr_valid(self.pgp) then begin
    n=n_elements(*self.pgp)
    result=make_array(n)
    for i=0,n-1 do result[i]=(*self.pgp)[i]->v0()
    return, result
  end else return, -1
end

;normally locates all points within alog10(delta_t)<.1 and delta_v<10 km/s (this should be only one)
;/next_point forces to locate the closest point already existing in the grid, whatever the differences in v0,n0 are
function model_grid::findgp,n0,v0,next_point=next_point
  if v0 lt 1e6 then v0*=1e5
  if n0 lt 1e3 then n0=10.^n0
  if keyword_set(next_point) then begin
    result=where(abs(alog10(self->n0())-alog10(n0)) eq min(abs(alog10(self->n0())-alog10(n0))) and abs(self->v0()-v0) eq min(abs(self->v0()-v0))) 
   endif else begin
    result=where(abs(alog10(self->n0())-alog10(n0)) lt .1 and abs(self->v0()-v0) lt 1e6)
  endelse
  if n_elements(result) ne 1 then message,'undefined result, multiple mathces to n0,v0?'+string(n0)+', '+sting(v0)
  return, result
end

function model_grid::addgp,n0,v0
  if ptr_valid(pgp) && (self->findgp(n0,v0))[0] ne -1 then begin
     message,'gridpoint with params (n0,v0) already exists'+string(n0)+' '+string(v0)
  endif else begin
    newgp=obj_new('grid_point',n0,v0)
    self.pgp=ptr_valid(self.pgp) ? ptr_new([*self.pgp,newgp]) : ptr_new(newgp)
  endelse
  return,self->findgp(n0,v0)
end

pro model_grid::readanalysegrid,filename
  restore, filename
  for i=0,n_elements(data.n0)-1 do begin
    match=self->findgp(data.n0[i],data.v0[i])
    ;this could later be extended to add data to already existing points, but for the time beeing this is more like an initialisation of a model_grid with the precomputed data from the old routines, so it shouldfind an empty field.
    if match[0] eq -1 then begin
      indgp=self->addgp(data.n0[i], data.v0[i])
      (*self.pgp)[indgp]->sethydparam,data.Te_0[i],data.Tion_0[i],data.maxdepth[i],data.Te_max[i],data.Tion_max[i], data.n_steps[i], data.diff2eq[i]
      if tag_exist(data,'o7r') then (*self.pgp)[indgp]->newline,8,7,21.60,data.o7r[i]
      if tag_exist(data,'o7i') then (*self.pgp)[indgp]->newline,8,7,21.80,data.o7i[i]
      if tag_exist(data,'o7f') then (*self.pgp)[indgp]->newline,8,7,22.10,data.o7f[i]
      if tag_exist(data,'ne9r') then (*self.pgp)[indgp]->newline,10,9,13.46,data.ne9r[i]
      if tag_exist(data,'ne9i') then (*self.pgp)[indgp]->newline,10,9,13.56,data.ne9i[i]
      if tag_exist(data,'ne9f') then (*self.pgp)[indgp]->newline,10,9,13.70,data.ne9f[i]
      if tag_exist(data,'olya') then (*self.pgp)[indgp]->newline,8,8,18.97,data.olya[i]
      if tag_exist(data,'olyb') then (*self.pgp)[indgp]->newline,8,8,16.01,data.olyb[i]
      if tag_exist(data,'nelya') then (*self.pgp)[indgp]->newline,10,10,12.13,data.nelya[i]
      if tag_exist(data,'nlya') then (*self.pgp)[indgp]->newline,7,7,24.78,data.nlya[i]
      if tag_exist(data,'n6r') then (*self.pgp)[indgp]->newline,7,6,28.79,data.n6r[i]
      if tag_exist(data,'n6i') then (*self.pgp)[indgp]->newline,7,6,29.10,data.n6i[i]
      if tag_exist(data,'n6f') then (*self.pgp)[indgp]->newline,7,6,29.54,data.n6f[i]
      if tag_exist(data,'fe17l1501') then (*self.pgp)[indgp]->newline,26,17,15.01,data.fe17l1501[i]
      if tag_exist(data,'fe17l1705') then (*self.pgp)[indgp]->newline,26,17,17.05,data.fe17l1705[i]
      if tag_exist(data,'fe17l1709') then (*self.pgp)[indgp]->newline,26,17,17.09,data.fe17l1709[i]
      if tag_exist(data,'fe17l1678') then (*self.pgp)[indgp]->newline,26,17,16.78,data.fe17l1678[i]
      if tag_exist(data,'fe17l1526') then (*self.pgp)[indgp]->newline,26,17,15.26,data.fe17l1526[i]
      if tag_exist(data,'fe18l1420') then (*self.pgp)[indgp]->newline,26,18,14.20,data.fe18l1420[i]
      if tag_exist(data,'fe18l1600') then (*self.pgp)[indgp]->newline,26,18,16.00,data.fe18l1600[i]
      if tag_exist(data,'fe18l1607') then (*self.pgp)[indgp]->newline,26,18,16.07,data.fe18l1607[i]
      if tag_exist(data,'fe18l1467') then (*self.pgp)[indgp]->newline,26,18,14.67,data.fe18l1467[i]
      if tag_exist(data,'fe19l1352') then (*self.pgp)[indgp]->newline,26,19,13.52,data.fe19l1352[i]
      if tag_exist(data,'fe19l1380') then (*self.pgp)[indgp]->newline,26,19,13.80,data.fe19l1380[i]
      if tag_exist(data,'fe19l1466') then (*self.pgp)[indgp]->newline,26,19,14.66,data.fe19l1466[i]
    endif else message,'gridpoint with params (n0,v0) already exists'+string(data.n0[i])+' '+string(data.v0[i])
  endfor
end

function model_grid::getlineflux, iz,ion,lambda,dlambda,n0,v0
case n_params() of
  4: begin
    result=make_array(n_elements(*self.pgp))
    for i=0,n_elements(*self.pgp)-1 do result[i]=(*self.pgp)[i]->getlineflux(iz,ion,lambda,dlambda)
    return, result
   end
  5: begin
    index=n0
    if index le n_elements(self->n0()) then begin
      return,(*self.pgp)[index]->getlineflux(iz,ion,lambda,dlambda)
    endif else begin
      message,'Index out of range'
      return,-1
    endelse
   end
  6: begin
   index=self->findgp(n0,v0)
   if index[0] eq -1 then begin
     message,'No n0,v0 gridpoint with values'+string(n0)+' , '+string(v0)
     return,-1
   endif else return, (*self.pgp)[index]->getlineflux(iz,ion,lambda,dlambda)
   end
  else: begin
    Message,'calling sequence model_grid->getline(iz,ion,lambda,dlambda[,n0,v0]) or ...,lambda,dlambda[,index])'
    return, -1
   end
endcase
end

function model_grid::INIT,abund_file
  common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
  read_abund, abund_file,abund,abund_ref
  return,1
end

pro model_grid::CLEANUP
  if ptr_valid(self.pgp) then for i=0,n_elements(*self.pgp)-1 do (*self.pgp)[i]->cleanup
  ptr_free,self.pgp
end

pro model_grid__define
  void={model_grid,pgp:ptr_new()} ;pgp: Pointer to Grid Points
end

