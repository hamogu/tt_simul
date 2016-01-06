pro tt_simul_gridnew, starti,startj=startj
common tt_units
restore,'/data/hspc62/steh305/idl/tt_results/grid/rad6e3newfe.sav'

if n_elements(startj) eq 0 then startj=0
v=(findgen(17)*25.+200.)*km/s
n0=10.^(findgen(13)/3.+10.)* 2.08316e-24
for i=starti,16 do begin
  for j=startj,12 do begin  
    test_index=where(data.v0 eq v[i] and abs(alog10(data.n0*2.08316e-24)-alog10(n0[j])) lt 0.1)
    if ((data.maxdepth)[test_index] gt 800*1e5) then begin
      print, 'Working on number ',i,j
      print, 'v0=',v[i],'n0=',n0[j]
      tt_simul,v[i],n0[j],2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/grid/deep/'+strcompress(string(i)+'_'+string(j),/remove_all),/quiet
    endif  
  endfor
  startj=0
endfor  
end