;this procedure automatically runs a grid of models
;it needs the following:
;in the working directory to_calculate.sav with to_calculate=int[17,13] eq 1 where simulations shall be performed, 0 otherwise
;in the working directory the abund file to be used , called 'tt_simul.abund'


pro tt_simul_list
common tt_units
@~/idl_startup.pro ;sets path to Chianti database
!abund_file ='tt_simul.abund'
!except=0
v=(findgen(40)*25.+200.)*km/s
n0=10.^(findgen(22)/3.+7.)* 2.08316e-24
n0_index=indgen(22)-9


repeat begin
for i=0,39 do begin
  for j=0,21 do begin  
      ;test if calculation is necessary      
      restore, 'to_calculate.sav'
      if to_calculate[i,j] eq 1 then begin
        if ~file_test(strcompress(string(i)+'_'+string(n0_index[j]),/remove_all)+'.sav') then begin ;test file for existence
          fileinfo=file_info(strcompress(string(i)+'_'+string(n0_index[j]),/remove_all)+'.log')
	  if systime(1)-fileinfo.ctime gt 3e4 then begin
            print, 'Working on number ',i,j
            print, 'v0=',v[i],'n0=',n0[j]
            tt_simul,v[i],n0[j],2e4,2e4,'/data/hspc62/steh305/idl/tt_data/arnaud_rothenflug_ext43.ioneq',+strcompress(string(i)+'_'+string(n0_index[j]),/remove_all),/quiet,t_star=.1,max_number_of_steps=1e5,max_depth=1e20
	    restore, 'to_calculate.sav'
	    to_calculate[i,j]=0 
            save,to_calculate,file='to_calculate.sav'
          endif ;calculation
        endif else begin
          print, '.sav files exists, but to_calculate was not reset for ',i,j
	  restore, 'to_calculate.sav'
	  to_calculate[i,j]=0 
          save,to_calculate,file='to_calculate.sav'
        endelse
      endif
  endfor 
endfor
wait, 180
restore, 'to_calculate.sav'  
endrep until total(to_calculate) eq 0
end
