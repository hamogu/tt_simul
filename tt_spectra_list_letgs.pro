;somewhere in CHIANTI there is an "execute" which I need to find.

pro tt_spectra_list_letgs
common tt_units
print, 'aaaaa'
@~/idl_startup.pro
print,'bbb'
abundfile ='tt_simul.abund'
!abund_file = abundfile
!except=0
read_abund, abundfile, abund,abund_ref

if n_elements(startj) eq 0 then startj=4
v=(findgen(17)*25.+200.)*km/s
n0=10.^(findgen(13)/3.+10.);* 2.08316e-24

sngl_ion={master:'/data/hspc62/steh305/idl/tt_results/twhya/master2.ions',C:['C_6','C_5','C_4','C_4d','C_5d','C_3','C_2'],N:['N_6','N_7','N_4','N_3','N_2'],O:['O_7','O_8','O_6','O_6d','O_7d','O_5','O_4','O_3'],Neon:['Ne_9','Ne_10','Ne_8','Ne_8d','Ne_9d','Ne_7','Ne_6','Ne_5'],Mg:['Mg_11','Mg_12','Mg_10','Mg_10d','Mg_9','Mg_8','Mg_7','Mg_6'],Si:['Si_13','Si_14','Si_12','Si_12d','Si_11','Si_10','Si_9','Si_8','Si_7'],S:['S_10','S_11','S_12','S_13','S_14','S_14d','S_9','S_8','S_7','S_6','S_5'],Fe:['Fe_17','Fe_18','Fe_19','Fe_20','Fe_21','Fe_22','Fe_23','Fe_24','Fe_17d','Fe_18d','Fe_19d','Fe_20d','Fe_21d','Fe_22d','Fe_23d','Fe_24d','Fe_16','Fe_15','Fe_14','Fe_13','Fe_12','Fe_11','Fe_16d','Fe_15d','Fe_14d','Fe_13d','Fe_12d','Fe_11d']}
naddparm=8

min_lambda=6.
max_lambda=170.
ang2kev=12.39854

lambda_scale=findgen(32800)*.005+6.
n_energies=32800

base_abund='/usr/local/hssoft/chianti/dbase/abundance/version_3/grevesse_sauval98.abund'

print, 'tt_spectra_list_letgs is starting'

repeat begin
for i=0,16 do begin
  for j=0,12 do begin  
      ;test if calculation is necessary      
      restore, 'to_calculate.sav'
      if to_calculate[i,j] eq 1 then begin
        if ~file_test(strcompress(string(i)+'_'+string(j),/remove_all)+'spectrum_letgs.sav') then begin
          fileinfo=file_info(strcompress(string(i)+'_'+string(j),/remove_all)+'spectrum_letgs.sav')
	  if systime(1)-fileinfo.ctime gt 3e4 then begin
	    note='Working on this case since'+systime()
            save,note,file=strcompress(string(i)+'_'+string(j),/remove_all)+'spectrum_letgs.sav'
            print, 'Working on number ',i,j
            print, 'v=',v[i],'n0=',n0[j]
            base_spectra={PARAMVAL:[0.,0.],INTPSPEC:fltarr(n_energies)}
	    for k=1,naddparm do base_spectra=add_tag(base_spectra,fltarr(n_energies),'ADDSP00'+strtrim(string(k),2))
   	    base_spectra.paramval=[v[i],n0[j]]

	    restore,strcompress(string(i)+'_'+string(j),/remove_all)+'.sav'
    	    prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens
    	    elecdens=reform(elecdens)
	    ionfile=strcompress(string(i)+'_'+string(j),/remove_all)+'.ioneq'
    	    read_ioneq,ionfile,ioneq_logt,ioneq,ioneq_ref
    	    logem=reform(alog10(volem))
    
  	    tt_synthetic, min_lambda,max_lambda,/tt_s,logem=logem,logt_i=ioneq_logt,dens=elecdens, /phot,rphot=1,radtem=6e3, ioneq_n=ionfile,output=lines, master=sngl_ion.master
	    make_chianti_spec, lines,lambda_scale,spectrum,/cont,abund=abundfile,/phot
            base_spectra.intpspec=spectrum.spectrum;*factor
    	    for elem=1,naddparm do begin
    	      tt_synthetic, min_lambda,max_lambda,/tt_s,logem=logem,logt_i=ioneq_logt, dens=elecdens,/phot,rphot=1,radtem=6e3,ioneq_n=ionfile,output=lines,sngl=sngl_ion.(elem)
      	      make_chianti_spec, lines,lambda_scale,spectrum,abund=base_abund,/phot
      	      base_spectra.(elem+1)=spectrum.spectrum;*factor
    	    endfor
    	    save,base_spectra,i,j,file=strcompress(string(i)+'_'+string(j),/remove_all)+'spectrum_letgs.sav'
            restore,'to_calculate.sav'
            to_calculate[i,j]=0 ;test file for existence
            save,to_calculate,file='to_calculate.sav'
	  endif
        endif
      endif
  endfor 
endfor
wait, 180
restore, 'to_calculate.sav'  
endrep until total(to_calculate) eq 0
end

