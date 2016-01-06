function findmax, array
maxpos=make_array(n_elements(array[*,0]))
for i=0,n_elements(array[*,0])-1 do begin
  temp=max(array[i,*],maxpo)
  maxpos[i]=maxpo
endfor  
return, maxpos
END
