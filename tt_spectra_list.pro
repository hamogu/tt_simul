pro tt_spectra_list_write_elem, elements,abundfile,base_abund
element=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si', $
         'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co', $
         'Ni','Cu','Zn']
read_abund,abundfile,abund,abund_ref
; --- Write new abund file ---
fileinfo=file_info('temp_master.abund')
if ~(fileinfo.exists) then begin
  openw,af,'temp_master.abund',/get_lun
  for i=0,29 do if abund[i] gt 0. and where(strupcase(elements) EQ strupcase(element[i])) eq -1 then printf,af,i+1,alog10(abund[i])+12.,element[i], format='(i4,f8.2,a6)'
  printf,af,'-1'
  printf,af,'%file:' ,'temp_master.abund'
  printf,af,'%intermedient file for generation of single ion continua; Moritz Guenther'
  printf,af,'% file can be deleted'
  printf,af,'-1'
  free_lun,af
endif

read_abund,base_abund,abund,abund_ref
for i=0,n_elements(elments)-1 do begin
  iz=where(strupcase(elements[k]) EQ strupcase(element))
  ; --- Write new abund file ---
  fileinfo=file_info('temp_'+strtrim(elements[i],2)+'.abund')
  if ~(fileinfo.exists) then begin
    openw,af,'temp_'+strtrim(elements[i],2)+'.abund',/get_lun
    printf,af,iz+1,alog10(abund[iz])+12.,element[iz], format='(i4,f8.2,a6)'
    printf,af,'-1'
    printf,af,'%file:' ,'temp_elem.abund'
    printf,af,'%intermediate file for generation of single ion continua; Moritz Guenther'
    printf,af,'% file can be deleted'
    printf,af,'-1'
    free_lun,af
  endif
endfor

end

pro tt_spectra_list_clean_elem, elements
fileinfo=file_info('temp_master.abund')
if (fileinfo.exists) then file_delete,'temp_master.abund'
for i=0,n_elements(elments)-1 do begin
  fileinfo=file_info('temp_'+strtrim(elements[i],2)+'.abund')
  if (fileinfo.exists) then file_delete,'temp_'+strtrim(elements[i],2)+'.abund'
endfor
end


pro tt_spectra_list
common tt_units
@~/idl_startup.pro

radtemp=1.

abundfile ='tt_simul.abund'
!abund_file = abundfile
!except=0
;if n_elements(startj) eq 0 then startj=4
v=(findgen(40)*25.+200.)*km/s
n0_index=indgen(22)-9
n0=10.^(findgen(22)/3.+7.);* 2.08316e-24


sngl_ion={master:'/data/hspc62/steh305/idl/tt_results/twhya/master2.ions',C:['C_6','C_5','C_4','C_4d','C_5d','C_3','C_2'],N:['N_6','N_7','N_4','N_3','N_2'],O:['O_7','O_8','O_6','O_6d','O_7d','O_5','O_4','O_3'],Neon:['Ne_9','Ne_10','Ne_8','Ne_8d','Ne_9d','Ne_7','Ne_6','Ne_5'],Mg:['Mg_11','Mg_12','Mg_10','Mg_10d','Mg_9','Mg_8','Mg_7','Mg_6'],Si:['Si_13','Si_14','Si_12','Si_12d','Si_11','Si_10','Si_9','Si_8','Si_7'],S:['S_10','S_11','S_12','S_13','S_14','S_14d','S_9','S_8','S_7','S_6','S_5'],Fe:['Fe_17','Fe_18','Fe_19','Fe_20','Fe_21','Fe_22','Fe_23','Fe_24','Fe_17d','Fe_18d','Fe_19d','Fe_20d','Fe_21d','Fe_22d','Fe_23d','Fe_24d','Fe_16','Fe_15','Fe_14','Fe_13','Fe_12','Fe_11','Fe_16d','Fe_15d','Fe_14d','Fe_13d','Fe_12d','Fe_11d']}
naddparm=8
elements=['C','N','O','Ne','Mg','Si','S','Fe']

min_lambda=1.
max_lambda=50.
ang2kev=12.39854
kev_scale=findgen(12000)*0.00025+0.000125
lambda_scale=ang2kev/kev_scale
lambda_scale=lambda_scale[where(lambda_scale ge min_lambda and lambda_scale le max_lambda)]
n_energies=11008

base_abund=!xuvtop +'/abundance/version_3/grevesse_sauval98.abund'


tt_spectra_list_write_elem, elements,abundfile,base_abund


repeat begin
for i=0,39 do begin
  for j=0,21 do begin  
      ;test if calculation is necessary      
      restore, 'to_calculate_spectra.sav'
      if to_calculate[i,j] eq 1 then begin
        restore,'to_calculate.sav'
        if to_calculate[i,j] eq 0 then begin ;in case this runs simulatneously with the tt_simul, check that the i_j.sav is already finished
        ;if file_test(strcompress(string(i)+'_'+string(j),/remove_all)+'.sav') then begin  
          fileinfo=file_info(strcompress(string(i)+'_'+string(n0_index[j]),/remove_all)+'spectrum.sav')
	  if systime(1)-fileinfo.ctime gt 3e4 then begin
	    note='Working on this case since'+systime()
            save,note,file=strcompress(string(i)+'_'+string(n0_index[j]),/remove_all)+'spectrum.sav'
            print, 'Working on number ',i,j
            print, 'v=',v[i],'n0=',n0[j]
            base_spectra={PARAMVAL:[0.,0.],INTPSPEC:fltarr(n_energies)}
	    for k=1,naddparm do base_spectra=add_tag(base_spectra,fltarr(n_energies),'ADDSP00'+strtrim(string(k),2))
   	    base_spectra.paramval=[v[i],n0[j]]

	    restore,strcompress(string(i)+'_'+string(n0_index[j]),/remove_all)+'.sav'
    	    prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens
    	    elecdens=reform(elecdens)
	    ionfile=strcompress(string(i)+'_'+string(n0_index[j]),/remove_all)+'.ioneq'
    	    read_ioneq,ionfile,ioneq_logt,ioneq,ioneq_ref
    	    logem=reform(alog10(volem))
    
  	    tt_synthetic, min_lambda,max_lambda,/tt_s,logem=logem,logt_i=ioneq_logt,dens=elecdens, /phot,rphot=1,radtem=radtemp, ioneq_n=ionfile,output=lines, master=sngl_ion.master
	    make_chianti_spec, lines,lambda_scale,spectrum,/cont,abund=abundfile,/phot
	    
            base_spectra.intpspec=spectrum.spectrum;*factor
    	    for elem=1,naddparm do begin
    	      tt_synthetic, min_lambda,max_lambda,/tt_s,logem=logem,logt_i=ioneq_logt, dens=elecdens,/phot,rphot=1,radtem=radtemp,ioneq_n=ionfile,output=lines,sngl=sngl_ion.(elem)
      	      make_chianti_spec, lines,lambda_scale,spectrum,abund='temp_'+strtrim(elements[i],2)+'.abund',/phot, /cont
	      ;Hydrogen cannot be switched of, because the read_abund scales everything realtiv to hydrogen
	      ;here its abund is noted with 0 instead of 12, so the abund of the elements in question relativ to that is 10^12 higher than normal
	      ;therefore the factor 1e-12, in total this suppresses the H continuum by 12 orders of magnitude
      	      base_spectra.(elem+1)=spectrum.spectrum*1e-12;*factor
    	    endfor
    	    save,base_spectra,i,j,file=strcompress(string(i)+'_'+string(n0_index[j]),/remove_all)+'spectrum.sav'
            restore,'to_calculate_spectra.sav'
            to_calculate[i,j]=0 
            save,to_calculate,file='to_calculate_spectra.sav'
	  endif
        endif
      endif
  endfor 
endfor
wait, 180
restore, 'to_calculate_spectra.sav'  
endrep until total(to_calculate) eq 0
tt_spectra_list_clean_elem, elements
end

