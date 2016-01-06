function v0n02gridfile,v0,n0
if v0 gt 1e6 then v0=v0/1e5
if n0 gt 1e3 then n0=alog10(n0)
v=findgen(30)*25+200
n=findgen(20)/3.+10
testv=min(abs(v-v0),v_ind)
testn=min(abs(n-n0),n_ind)
if testv gt 10. then message,'Can not find matching v0'
if testn gt .1 then message,'Can not find matching n0'
return,strtrim(string(v_ind),2)+'_'+strtrim(string(n_ind),2)
end
