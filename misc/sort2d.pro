;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	SORT"D
;
; PURPOSE:
;	This function transforms one vector into an 2-dim array. This is meant to be data of the form z=z(x,y). Two input vectors of the same size provide the x and y.
;
; CATEGORY:
;	misc
;
; CALLING SEQUENCE:
;	Result = SORT2D(vector,index1,index2)
;
; INPUTS:
;	vector:		The vector to be sorted
;	index1:		A vector of the same size, it provides the X values for sorting
;	index2:		A vector of the same size, it provides the Y values for sorting
;;
; EXAMPLE:
;		F = SORT2D([1,2,3,4],[1,1,2,2],[1,2,1,2])
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Günther, 28.9.2005
;-


function sort2d, array, index1,index2

n1=n_elements(uniq(index2,sort(index2)))
n2=n_elements(uniq(index1,sort(index1)))

if n1*n2 ne n_elements(array) then begin
  print, 'index vectors do not contain the correct number of unique entires to generate an 2d array correctly'
  return, -1
endif
  
i1=sort(index2)
index22d=index1[i1]
array2d=array[i1]

array2d=reform(array2d,[n2,n1],/overwrite)
index22d=reform(index22d,[n2,n1],/overwrite)

for i=0,n1-1 do begin
   i2=sort(index22d[*,i])
   array2d[*,i]=array2d[i2,i]
endfor  
return, array2d
end