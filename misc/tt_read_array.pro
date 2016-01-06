;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_READ_ARRAY
;
; PURPOSE:
;	This procedure reads a generic float array from a file
;	The first line must give the dimensions of the array e.g. 2 4 4
;	The rest contains the data in the form index1 index2 data1 data2 data3 data4
;	e.g. 1 2 3. 4. 5. 6. corresponds to arr[1,2,*]=[3.,4.,5.,6.]
;	works for any number of dimensions
;	references or additional information can be given at the bottom
;	enclosed in two lines "-1"
;
; CAVEATS:
;	IDL array idices start with 0, if the data starts with 1 specify all dimensions except the 
;	last one higher, so read 3*4*5 data in a 4*5*5 array
;	BE CAREFUL: this means data point [1,1,1] is in coeff[1,1,0] !
;
; CATEGORY:
;	files
;
; CALLING SEQUENCE:
;	tt_read_array ,filename,coeff,ref
;
; INPUTS:
;	filename:	will be selected with a widget if not specified
;
; OUTPUTS:
;	coeff:		the array
;	ref:		an array of references
;
; TEST STATUS:
;	tested with real data of 3 dimensions
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther 31.01.2005
;-


pro tt_read_array ,filename,coeff,ref
;
;
if n_params(0) lt 3 then begin
   print,''
   print,' >   tt_read_array ,filename,arr,ref'
   print,'      or'
   print,' > read_ioneq,''',''', arr, ref'
   print,''
   return
endif
;
if strtrim(filename,2) eq '' then begin
;
    dir=concat_dir(!xuvtop,'ioneq')
    filename=dialog_pickfile(path=dir,filter='*.dat',title='Select Array File')
    print,' selected:  ',filename
endif
;
;
;
;
openr,lu,filename,/get_lun
;
string1=' '
str=''
;
;
readf,lu,str  ; read dimensions
str=strtrim(str,2)
dimensions=strsplit(str,/extract)
n_dim=n_elements(dimensions)
coeff=fltarr(dimensions)
formatstr='('+strtrim(n_dim-1,2)+'i3,'+strtrim(dimensions[n_dim-1],2)+'e10.2)'
indexline=intarr(n_dim-1)
dataline=fltarr(dimensions[n_dim-1])

;
;
while strpos(string1,'-1') EQ  -1 or strpos(string1,'-1') GT 2  do begin
readf,lu,string1
if(strpos(string1,'-1')   EQ  -1 or strpos(string1,'-1') GT 2) then begin

  reads,string1,indexline, dataline			;,format=formatstr
  ;coeff[indexline,*]=dataline[*]
  if n_dim eq 2 then coeff[indexline[0],*]=dataline[*]
  if n_dim eq 3 then coeff[indexline[0],indexline[1],*]=dataline[*]
  if n_dim eq 4 then coeff[indexline[0],indexline[1],indexline[2],*]=dataline[*]
  if n_dim eq 5 then coeff[indexline[0],indexline[1],indexline[2],indexline[3],*]=dataline[*]
  if n_dim eq 6 then coeff[indexline[0],indexline[1],indexline[2],indexline[3],indexline[4],*]=dataline[*]
  if n_dim eq 7 then coeff[indexline[0],indexline[1],indexline[2],indexline[3],indexline[4],indexline[5],*]=dataline[*]
  if n_dim eq 8 then coeff[indexline[0],indexline[1],indexline[2],indexline[3],indexline[4],indexline[5],indexline[6],*]=dataline[*]

endif
endwhile
;
;  get references
refstring=strarr(100)
nref=0
;
string1=' '
while strpos(string1,'-1') EQ  -1  do begin
readf,lu,string1


if(strpos(string1,'-1') EQ -1) and (strpos(string1,'%file') EQ -1) then begin
  refstring(nref)=string1
  nref=nref+1
endif
endwhile
;
ref=refstring(0:nref-1)
;
;
free_lun,lu
;
;
end
