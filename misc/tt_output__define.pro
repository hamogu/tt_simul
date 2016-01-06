;+
; PROJECT:
;	TT_SIMUL
;
; CLASS_NAME:
;	tt_output
;
; PURPOSE:
;	A TT_OUTPUT object organizes the results of one tt_simul run. Calculated iononisation values and hydrodynamic 
;	variables are stored temporarily and can be written to the disk later. THe ioneq data is stored in CHIANTI *.ioneq format
;	to allow further analysis with the CHIANTI package.
;
; CATEGORY:
;	data management
;
;
; CREATION:
;       On creation an TT_OUTPUT object is initialized with the path to store the data including the filemane without extension,
;	optionally the current program version and the abundance fileif this is not stored in the system varable !abund_file 
;	Example: obj2=obj_new('tt_output','./test','0.44') 
;
; METHODS:
;	TT_Output objects have internal methods, but IDL does not provide a mechanism for method encapsulation. Internal methods
;	are not listed here. They are documented only directly in the code.
;       tt_output::write_data
;	tt_output::add_data, ioneq,T,hydrodyn
;	tt_output::get_hydrodyn,hydrodyn
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther, 21.03.2005

;
; =============================================================
; Description of methods
; =============================================================
;
; METHODNAME:
;       tt_output::add_data
;
; PURPOSE:
;       Use this method to add data.
;
; CALLING SEQUENCE:
;	Obj -> add_data, ioneq,T,hydrodyn
;
; INPUTS:
;	ioneq:	an 2-D array containig the ionisation state. THis should be a 2-dim version of the ioneq
;		from the CHINATI common block common elements
;		Note: refine this to a 2-D array for a single temperature only, e.g add_data, reform(ioneq[0,*,*])
;	T:	temperatue in K. Saved in the .ioneq file, so will be used by CHIANTI as temperature for spectrum
;		calculations
;	hydrodyn: An array containing hydrodynamical data, the format is not explicitly specified and my vary with program version.
;		e.g. [t_ion, T_e,velocity,...]
;
; EXAMPLE:
;	Add some simple (non realistic) data.
;	obj->add_data,[[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]],10000.,[21,34]
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther, 21.03.2005
;
; METHODNAME:
;       tt_output::write_data
;
; PURPOSE:
;       Use this method to write data to the disk. Several files are produced. For the important ones (*.hydrodyn and *.ioneq)
;	it is checke if they exist and they are backuped and overwritten, if the backup already exits it is deleted.
;	These files serve two purposes:
;	1) keep data for analysis in general (*.hydrodyn,*.ioneq,*.master - IDL format)
;	2) divide the simulation in two regions where the electron temperature is a monotone function
;	   This is meant for use with CHIANTI to get complete spectra (*.dem,*.ioneq,*.tene)
;
; CALLING SEQUENCE:
;	Obj -> write_data
;
; EXAMPLE:
;	write to disk
;	obj->write_data
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther, 21.03.2005
;
; METHODNAME:
;       tt_output::get_ionstate
;
; PURPOSE:
;       Use this method get the array with ioneq values for a single ion, so you can e.g. plot abundance of O VII to x. If no data exists the method
;	returns -1
;
; CALLING SEQUENCE:
;	Obj -> get_ionstate,iz,ion,ionsate
;
; INPUTS:
;	iz:	atomic number	e.g. 8 for oxygen
;	ion:	ion in spectroscopic notation	e.g. 7 for O VII
;
; OUTPUTS:
;	ionstate: An array containing ioneq datain the same order as in hydrodyn
;
; EXAMPLE:
;	get data for O VII
;	obj->get_ionstate, 8,7,ionstate
;
; METHODNAME:
;       tt_output::get_hydrodyn
;
; PURPOSE:
;       Use this method get the array with hydrodynamic data directly without writing to disk first. If no data exists the method
;	returns -1
;
; CALLING SEQUENCE:
;	Obj -> get_hydrodyn,hydrodyn
;
; OUTPUTS:
;	hydrodyn: An array containing hydrodynamical data, the format is not explicitly specified and may vary with program version.
;		e.g. [t_ion, T_e,velocity,...]
;
; EXAMPLE:
;	get hydorynamic data
;	obj->get_hydrodyn, hydrodyn
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther, 21.03.2005
;
; METHODNAME:
;       tt_output::log
;
; PURPOSE:
;       adds some text tot the logfile
;
; CALLING SEQUENCE:
;	Obj -> log,string
;
; INPUTS:
;	text: a string which is written in the .log file
;
; EXAMPLE:
;	obj->log,'Log this text'
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther, 29.03.2005
;
; METHODNAME:
;       tt_output::plog
;
; PURPOSE:
;       adds some text to the logfile and prints it to the screen
;
; CALLING SEQUENCE:
;	Obj -> plog,string
;
; INPUTS:
;	text: a string which is written in the .log file
;
; EXAMPLE:
;	obj->plog,'Log this text'
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther, 01.04.2005
;-


;ONLY for internal calls
pro tt_output::add_t,T
  *self.plogt= (n_elements(*self.plogt)eq 0) ? alog10(T): [alog10(T),*self.plogt]
end
;ONLY for internal calls
;Not all elemtents are actually used in most situations. This routine filters out the elements which do appear in the ioneq
;an array corresponding to the first and second column in the CHIANTI *.ioneq files is stored in *self.pioneqindex
pro tt_output::define_saveform,ioneq
  ;create an index of elements whih are used
  for i=0,n_elements(ioneq[*,0])-1 do if total(ioneq[i,*] gt 0) then index= (n_elements(index)eq 0) ? i+1: [index,i+1]
  ;build an array matching the first and second column in .ioneq files
  ;this cannot be done in an array operation since some elements from ioneq are jumped
  for i=0,n_elements(index)-1 do begin	;loop over elements
    for ion=1,index[i]+1 do begin	;loop over ions of that element ;+1 to get e.g. HI and HII
      *self.pioneqindex= (n_elements(*self.pioneqindex)eq 0) ? [index[i],ion]: [[*self.pioneqindex],[index[i],ion]]
    endfor  
  endfor  
  self.saveformdefined=1
end  
;ONLY for internal calls
;ioneq as 2-Dim array
pro tt_output::add_ioneq,ioneq
  if self.saveformdefined eq 0 then self->define_saveform,ioneq
  new_column=make_array(1,n_elements((*self.pioneqindex)[0,*]))
  ;
  for i=0,n_elements((*self.pioneqindex)[0,*])-1 do begin
    new_column[0,i]=ioneq[(*self.pioneqindex)[0,i]-1,(*self.pioneqindex)[1,i]-1]
  endfor
  *self.pioneq=(n_elements(*self.pioneq) eq 0) ? new_column : [new_column,*self.pioneq]
  self.saved=0
end
;for internal calls
pro tt_output::add_iondata, ioneq,T
  self->add_t,T
  self->add_ioneq,ioneq
  self.saved=0
end
;for internal calls
pro tt_output::add_hydrodyn, hydrodyn
  *self.phydrodyn=(n_elements(*self.phydrodyn) eq 0) ? [hydrodyn] : [[hydrodyn],[*self.phydrodyn]]
  self.saved=0
end  
;use this for adding data
pro tt_output::add_data, ioneq,T,hydrodyn
  self->add_hydrodyn,hydrodyn
  self->add_iondata,ioneq,T
  self.saved=0
end  
pro tt_output::get_hydrodyn,hydrodyn
  hydrodyn=(n_elements(*self.phydrodyn) eq 0)?  -1 : *self.phydrodyn
end
pro tt_output::get_ionstate,iz,ion,ionstate
  ionstate=-1
  if self.saveformdefined then begin
    row=where(((*self.pioneqindex)[0,*] eq iz) and ((*self.pioneqindex)[1,*] eq ion))
    if (row ne -1) then ionstate=(*self.pioneq)[*,row]
  endif  
end 
pro tt_output::write_comment,lu
  printf,lu,'%date:'+self.date
  printf,lu,'%program version:'+string(self.version)
  printf,lu,'%produced as output of TT-SIMUL, simulating the accretion shock of T Tauri Stars' 
  printf,lu,'%Hamburger Sternwarte, Moritz Guenther'
  printf,lu,'%data does not represent LTE'
end  
pro tt_output::write_ioneq_data,filename,ioneq,logt
  openw,lu,filename,/get_lun
  n_t=n_elements(logt)
  n_elem=max((*self.pioneqindex)[0,*])
  ;header
  printf, lu,format='(I3,I4)',n_t,n_elem
  ;in CHIANTI .ioneq data is sorted in ascending order of logt
  ;index=sort(*self.plogt)
  index=findgen(n_elements(logt))	;print unsortet .ioneq
  ;logt
  lineformat='('+string(n_t)+'(f8.4))'
  printf,lu,format=lineformat,logt[index]
  ;ioneqdata
  lineformat='(i3,i3,'+string(n_t)+'(e10.3))'
  for i=0,n_elements((*self.pioneqindex)[0,*])-1 do printf, lu,format=lineformat, (*self.pioneqindex)[*,i],ioneq[index,i]
  ;comments
  printf,lu,' -1'
  name=strsplit(filename,'/',count=i,/extract)
  printf,lu,'%filename:'+name[i-1]
  printf,lu,'%corresponding hydrodynfile:'+strmid(name[i-1],0,strlen(name[i-1])-6)+'.hydrodyn'
  self->write_comment,lu
  printf,lu,' -1'
  free_lun,lu
end
;for regular ioneq-files
pro tt_output::write_ioneq
  filename=self.path+'.ioneq'
  if file_test(filename)then begin
    self->log, filename+ "already exists. It is overwritten. Backup: .ioneq~"
    file_move,filename,filename+"~",/overwrite
  endif  
  ioneq=*self.pioneq
  logt=*self.plogt
  self->write_ioneq_data,filename,ioneq,logt
end
pro tt_output::write_hydrodyn
  filename=self.path+'.hydrodyn'
  if file_test(filename)then begin
    self->log, filename+"already exists. It is overwritten. Backup: .hydrodyn~"
    file_move,filename,filename+"~",/overwrite
  endif  
  openw,lu,filename,/get_lun
  ;data
  lineformat='('+string(n_elements((*self.phydrodyn)[*,0]))+'(e17.10))'
  for i=0,n_elements((*self.phydrodyn)[0,*])-1 do printf, lu,format=lineformat, (*self.phydrodyn)[*,i]
  ;comments
  printf,lu,' -1'
  name=strsplit(filename,'/',count=i, /extract)
  printf,lu,'%filename:'+name[i-1]
  printf,lu,'%corresponding ioneqfile:'+strmid(name[i-1],0,strlen(name[i-1])-9)+'.ioneq'
  self->write_comment,lu
  printf,lu,' -1'
  free_lun,lu
end
pro tt_output::write_master
  filename=self.path+'.sav'
  self->log, 'saving data with IDL save in'+filename
  ioneq=*self.pioneq
  ioneqindex=*self.pioneqindex
  hydrodyn=*self.phydrodyn
  save,ioneq,ioneqindex,hydrodyn,filename=filename
end
;for use with CHIANTI
pro tt_output::write_tene
  t_e=(*self.phydrodyn)[5,*]
  n_e=(*self.phydrodyn)[2,*]*(*self.phydrodyn)[6,*]
  maxtemp=max(t_e,index)
  ;write data past t_e maximum
  filename=self.path+'_2.tene'
  self->log, 'saving tene raw data in '+filename
  openw,lu,filename,/get_lun
  printf,lu,[t_e[0,0:index],n_e[0,0:index]]
  free_lun,lu
  ;write data between shock and T_e maximum
  filename=self.path+'_1.tene'
  self->log, 'saving tene raw data in '+filename
  openw,lu,filename,/get_lun
  ;the last 2 are marked with t=inf, because they represen input data
  printf,lu,[t_e[0,index:n_elements(t_e)-2],n_e[0,index:n_elements(t_e)-2]]
  free_lun,lu
end
;for use with CHIANTI
pro tt_output::write_ioneq_parts
  ioneq=*self.pioneq
  logt=*self.plogt
  t_e=(*self.phydrodyn)[5,*]
  maxtemp=max(t_e,index)
  ;write data past t_e maximum
  filename=self.path+'_2.ioneq'
  self->write_ioneq_data,filename,ioneq[0:index,*,*],logt[0:index]
  ;write data past t_e maximum
  ;the last 2 are marked with t=inf, because they represen input data
  filename=self.path+'_1.ioneq'
  self->write_ioneq_data,filename,ioneq[index:n_elements(t_e)-2,*,*],logt[index:n_elements(t_e)-2]
end
;for use with CHIANTI
pro tt_output::write_dem
  hydrodyn=*self.phydrodyn
  t_e=hydrodyn[5,*]
  n_steps=n_elements(hydrodyn[0,*])
  delta_te=-hydrodyn[5,*]+[[reform(max(hydrodyn[5,*]),1,1)],[hydrodyn[5,0:n_steps-2]]]
  delta_x=-hydrodyn[0,*]+[[reform(max(hydrodyn[0,*]),1,1)],[hydrodyn[0,0:n_steps-2]]]
  volem=delta_x*hydrodyn[2,*]^2*hydrodyn[6,*]
  dem=abs(volem/delta_te)
  maxtemp=max(t_e,index)
  ;write data past t_e maximum
  filename=self.path+'_2.dem'
  self->log, 'saving dem data past T_e maximum in '+filename
  openw,lu,filename,/get_lun
  printf,lu,[alog10(t_e[0,0:index]),alog10(dem[0,0:index])]
  printf,lu,' -1'
  self->write_comment,lu
  printf,lu,' -1'
  free_lun,lu
  ;write data between shock and T_e maximum
  filename=self.path+'_1.dem'
  self->log, 'saving dem data between shock and maximum in '+filename
  openw,lu,filename,/get_lun
  ;the last 2 are marked with t=inf, because they represen input data
  printf,lu,[alog10(t_e[0,index:n_elements(t_e)-2]),alog10(dem[0,index:n_elements(t_e)-2])]
  printf,lu,' -1'
  self->write_comment,lu
  printf,lu,' -1'
  free_lun,lu
end 
pro tt_output::write_data
  self->write_hydrodyn
  self->write_ioneq
  self->write_master
  self->write_tene
  self->write_dem
  self->write_ioneq_parts
  saved=1
end  
  
function tt_output::INIT, path, version, ABUND_FILE=abund_file
  case n_params() of
    1: self.version=''
    2: self.version=version
    else:Message,'Invalid number of Arguments'
  endcase     
  self.abundfile= (n_elements(abund_file) eq 0) ? !abund_file : abund_file
  self.date=systime()
  self.path=path
  self.pioneq=ptr_new(/Allocate_heap)
  self.pioneqindex=ptr_new(/Allocate_heap)
  self.plogt=ptr_new(/Allocate_heap)
  self.phydrodyn=ptr_new(/allocate_heap)
  return, 1
end
;checks if data was saved after last modification.
;if not is is written in default.ioneq and default.hydrodyn, an error is generated, but program execution continues
pro tt_output::CLEANUP
  if ~self.saved then begin
    self->log, "Destruction of Output Structure without saving"
    self->log, "Data is stored in current directory as default.*"
    self.path='./default'
    self->write_data
  endif  
  ;free heap space
  IF (PTR_VALID(self.pioneq)) THEN PTR_FREE,self.pioneq
  IF (PTR_VALID(self.pioneqindex)) THEN PTR_FREE,self.pioneqindex
  IF (PTR_VALID(self.plogt)) THEN PTR_FREE,self.plogt
  IF (PTR_VALID(self.phydrodyn)) THEN PTR_FREE,self.phydrodyn
end
pro tt_output::log, text
  filename=self.path+'.log'
  if ~file_test(filename) then begin
    openw, lu,filename,/get_lun
    printf,lu,'Log-file for tt_simul run'
    self->write_comment, lu
  endif else begin
    openw,lu,filename,/get_lun,/append
  endelse 
  printf,lu,text
  free_lun,lu
end
pro tt_output::plog, text
  print, text
  self->log,text
end  
pro tt_output__define
  void={tt_output,pioneq:ptr_new(),pioneqindex:ptr_new(),plogt:ptr_new(),phydrodyn:ptr_new(),version:'',date:'',path:'',abundfile:'',saved:0,saveformdefined:0}
end
