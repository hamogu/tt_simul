pro tt_simul_grid, starti,startj=startj
common tt_units

if n_elements(startj) eq 0 then startj=0
v=(findgen(17)*25.+200.)*km/s
n0=10.^(findgen(13)/3.+10.)* 2.08316e-24
for i=starti,16 do begin
  for j=startj,12 do begin
    print, 'Working on number ',i,j
    print, 'v0=',v[i],'n0=',n0[j]
    tt_simul,v[i],n0[j],2e4,2e4,'../tt_data/arnaud_rothenflug_ext43.ioneq','../tt_results/grid/'+strcompress(string(i)+'_'+string(j),/remove_all),/quiet
  endfor
  startj=0
endfor  
end