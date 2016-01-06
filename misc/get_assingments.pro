
;-------------- Parameter file handling  --------------------------------------------------------------

; get_assignments - Parsing of a parameter file for predefining variables (idl lacking proper constants).
;                   The par.file must use proper idl assignments, one per line, no if statements etc. allowed.
;
; The string array assigns is returned, it contains all succesfully parsed assignments.
; To perform the assignments on the mainprogram level, use eg.
;  assis=get_assignments("tst.dat", warn_level=0)
;  for tmpj=0,(n_elements(assis)-1) do  begin dummy=execute(assis[tmpj]) & print,format='(I4,$)',tmpj & endfor & print
;  delvar,assis    ;careful ! Seems to cause trouble.
; or more readable
;  for tmpj=0,(n_elements(assis)-1) do  begin
;    dummy=execute(assis[tmpj]) & print,format='(I4,$)',tmpj
;  endfor  &  print
;
;    file_str contains the filename of the parameter file, ="NONE" causes no parameter processing, ="" causes
;    automatically looking for a parameter file of the same name as the calling programme. ANYWAY the par-file
;    mustbe in the same directory as the calling program.
;  For warn_level not given or >9, a dummy assignment is performed to check the validity of the statements
;  (time consuming). An option would be handy, parsing the calling programme for assignments overriding
;  the defined 'constants' in the parameter file. Fatal conditions are always reported, even at warn_level=0.
;  Warn_level=-1 suppresseses everything but warnings.
;    maxassigns sets the max. number of assignments accepted.
;
function get_assignments, file_str, warn_level=warn_level, maxassings=maxassigns
if (keyword_set(maxassings) eq 0) then maxassigns=100                               ;set default max. no. of assignments
if (keyword_set(warn_level) eq 0) then warn_level=99                                ;set default warn-level
file_no=9
cm_char = ";"                                                   ;token for comments
mu_char = "&"                                                   ;token seperating multiple commands per line
routine_str ="get_assignments"
prgext_str = ".pro"                                  ;extension of program files
parext_str = ".par"                                  ;default extension to use for parameterfiles
tmpstruc = routine_info("$MAIN$", /source)  ;={name:"", path:""}

if (strpos(tmpstruc.path,"/") ne -1) then  pathsep_str = "/"  else pathsep_str = "\"       ;find out OS:
                                                                                           ;  path seperator LINUX vs.WIN

if (file_str eq "") then  begin                                        ;## file_str="" ? .. auto
  file_str = strmid( tmpstruc.path, $
                     0,(rstrpos(tmpstruc.path,prgext_str)-0) )  $  ;substitute prog.file extension
            + parext_str
endif else begin                                                       ;## no
  file_str = strmid( tmpstruc.path, $
                     0,(rstrpos(tmpstruc.path,pathsep_str)+1) )  $    ;get calling programs path incl. pathseparator
            + file_str
endelse

print
print, routine_str, file_str, $
       format='(">>> ",A ," processing... ", A)'
print, "Invoked by ", tmpstruc.path

if (file_str eq "NONE") then  $                           ;## file_str="NONE" ?
  return, cm_char+" No assignment file processed"        ;return just a comment for none

openr, file_no, file_str

assigns=make_array(maxassigns,/string)                                            ;create array

read_str =""                                                                      ;dummy assignment: set type to string
lineno=0  &  assino=0                                                             ;line counter, assignment counter
while (not eof(file_no)) do begin
  lineno = lineno+1
  readf, file_no, read_str
  read_str=strtrim(read_str ,2)                                                   ;->remove leading&trailing spaces
  if ((strpos(read_str,cm_char) ne 0) and (strlen(read_str) gt 2) ) then begin        ;# if line neither empty nor command
    before_comment=str_sep(read_str,cm_char)  &  before_comment=before_comment[0]       ;-> remove comments after command
  ;TBD: parse multiples here
    if ((strpos(read_str,"=") eq 0) or (strpos(read_str,mu_char) ne 0)) then begin       ;# if not assignment or multiple
  ;TBD: extract and keep name of variable assigned to
      if (assino ge maxassigns) then begin
          close, file_no                                                                      ;    close file
          print, lineno, format='("Too many assignments in line ", I3)'
          print, read_str
          message, "ERROR0 in "+routine_str                                              ;    report error
      endif
      assigns[assino] = before_comment
      assino=assino+1
      if (warn_level ge 0) then $
        print, assino,before_comment, format='(I3, ">", A, "<",$)'                     ;print assumed assignment
      if (warn_level ge 10) then begin                                               ;# dummy execute to check syntax ?
        ok_flag = execute(before_comment)
        if (ok_flag eq 0) then begin                                                     ;if not executeable...
          close, file_no                                                                      ;    close file
          print, lineno, format='("Non-executeable statement in line ", I3)'
          print, read_str
          message, "ERROR1 in "+routine_str                                              ;    report error
        endif
      endif
    endif else begin                                                                 ;# else report error
      close, file_no
      print, lineno, format='("Non-assignment or multiple statement in line ", I3)'
      print, read_str
      message, "ERROR2 in "+routine_str
    endelse
  endif
endwhile

print
print, routine_str, lineno, assino, $
       format='(">>> ",A ," processed ", I3, " lines, extracting ", I3, " assignments.")'
print
close, file_no

assigns=assigns[0:assino-1]                         ;reduce to nonempty
return, assigns
END