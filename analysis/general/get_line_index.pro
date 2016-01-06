;short documentaion
; This routine searches a structure for lines
; Ch_synthetic putputs a lines structure with different fields identifiing single lines
; for further processing it may be useful to find the index of lines
; snote : spectroscopic notation of ion "Ne X"
; wvl:    wavelength to be searched for [Angstroem]
; d_lambda: delta_lambda (true wavelength -- input wavelength), default=1.
; OUTPUT:
; -1:	ion not in databas
; -2: no line within d_lambda

function get_line_index, linearray,snote,wvl,d_lambda=d_lambda
if n_elements(d_lambda) ne 1 then d_lambda=1.
ion_index=where(linearray.snote eq snote)

if ion_index[0] eq -1 then begin  ;ion_index[0] not just ion_index, because cannot compare vector to -1
  return,-1
 endif else begin 
  found_min=where(abs(linearray[ion_index].wvl-wvl) le d_lambda)
  if found_min[0] eq -1 then return, -2 else return, ion_index[found_min]
endelse    
end