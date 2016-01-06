;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	CALCULATE_RATIOS
;
; PURPOSE:
;	This procedure calculates line ratios from the inpute data obs and simdata.
;
; CATEGORY:
;	fitting
;
; CALLING SEQUENCE:
;	CALCULATE_RATIOS, obs,simdata,obsratios,simratios,ratiofile=ratiofile
;
; INPUTS:
;	ratios:		a string array with the ratios to be computed, in a format compatible to calculate_ratios.pro
;			accepted format of a line is:
;			empty
;			;comments
;			(a+b)/c
;			a/b
;			a,b,c are the identifiers used as components of the data structure and as label in the observational data file
;	obs:		observatinal data from read_obs_data.pro
;	simdata:	output from analyse_grid.pro
;
; KEYWORDS:
;	ratiofile:	if set ratios is not a string array, but the path to a file with such an array written line by line
;
; OUTPUTS:
;	obsratios:	array with the calculates ratios
;			obsratios[*,0] contains the values, obsratio[*,1] the errors
;	simratios:	an array with the ratios for all model grid points
;
; TEST STATUS:
;	checked against by hand calculations
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 22.09.2005
;-

pro calculate_ratios, ratios,obs,simdata,obsratios,simratios,ratiofile=ratiofile

n_sim=n_elements(simdata.n0)

;read file which specifies the line ratios to be used
if keyword_set(ratiofile) then begin
  line='string'
  openr,lu,ratios,/get_lun
  while ~eof(lu) do begin
   readf,lu,line
   read_in=n_elements(read_in) gt 0 ? [read_in,line] : line
  endwhile
  free_lun, lu
  ratios=read_in
endif
  


obsratios=make_array(n_elements(ratios),2,val=0.)
simratios=make_array(n_elements(ratios),n_sim,val=0.)
for linenumber=0,(n_elements(ratios)-1) do begin
  
  case 1 of ;check if the string matches one of the supported formats - 1=true

    ;e.g (f+i)/r
    strmatch(ratios[linenumber],'(*+*)/*',/fold_case): begin
      ratio=strlowcase(strsplit(ratios[linenumber],'/',/extract))
      ;strmid cuts of the ( ) at beginning and end of string
      num=strlowcase(strsplit(strmid(ratio[0],1,strlen(ratio[0])-2),'+',/extract))
      numi1=where(strlowcase(obs.name) eq num[0])  ;index of first part in numarator
      numi2=where(strlowcase(obs.name) eq num[1])  ;index if second part in numarator
      dei=where(strlowcase(obs.name) eq ratio[1])  ;index of denominator
      ;check if all reqiered data exist in the data file
      if numi1 eq -1 or numi2 eq -1 or dei eq -1 then begin
       print, 'Ratio: ', ratios[linenumber],' needs data not contained in input file!'
       return
      endif  
      ;ratio
      obsratios[linenumber,0]=((obs.flux)[numi1]+(obs.flux)[numi2])/(obs.flux)[dei]
      ;error
      obsratios[linenumber,1]=sqrt(((obs.error)[numi1]^2.+(obs.error)[numi2]^2.)/(obs.flux)[dei]^2.+(((obs.flux)[numi1]+(obs.flux)[numi2])/(obs.flux)[dei]^2.)^2.*(obs.error)[dei]^2.)
      
      test=execute('simratios[linenumber,*]=(simdata.'+num[0]+'+simdata.'+num[1]+')/simdata.'+ratio[1])
      end
    ;e.g f/i
    strmatch(ratios[linenumber],'*/*',/fold_case): begin
      ratio=strlowcase(strsplit(ratios[linenumber],'/',/extract))
      numi=where(strlowcase(obs.name) eq ratio[0])    ;index of first part in numarator
      dei=where(strlowcase(obs.name) eq ratio[1])     ;index of denominator
      ;check if all reqiered data exist in the data file
      if numi eq -1 or dei eq -1 then begin
       print, 'Ratio: ', ratios[linenumber],' needs data not contained in input file!'
       return
      endif  
      ;ratio
      obsratios[linenumber,0]=(obs.flux)[numi]/(obs.flux)[dei]
      ;error
      obsratios[linenumber,1]=sqrt((obs.error)[numi]^2./(obs.flux)[dei]^2.+((obs.flux)[numi]/(obs.flux)[dei]^2.)^2.*(obs.error)[dei]^2.)
      
      test=execute('simratios[linenumber,*]=simdata.'+ratio[0]+'/simdata.'+ratio[1])
      end
    ;empty line
    ratios[linenumber] eq '': ;do nothing
    ;comment line
    strmatch(ratios[linenumber],';*'): ;do nothing
     
    else: print, 'Format of line ratio not supported: ',ratios[linenumber] ,'  --- line ignored' ;ignore line
  endcase 
endfor

end