pro tt_simul_newabund, starti,startj=startj
common tt_units
;!abund_file ='/data/hspc62/steh305/idl/tt_data/bptau.abund'
!abund_file ='/data/hspc62/steh305/idl/tt_data/twhya2.abund'
!except=0
if n_elements(startj) eq 0 then startj=4
v=(findgen(17)*25.+200.)*km/s
n0=10.^(findgen(13)/3.+10.)* 2.08316e-24
;v=[350.,450.,550.]*km/s
;n0=10.^[11.5,12.,12.5]* 2.08316e-24
for i=starti,15 do begin
  for j=startj,8 do begin  
    
      print, 'Working on number ',i,j
      print, 'v0=',v[i],'n0=',n0[j]
      tt_simul,v[i],n0[j],2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/twhyaabund2/'+strcompress(string(i)+'_'+string(j),/remove_all),/quiet
    
  endfor
  startj=4 
endfor  
end