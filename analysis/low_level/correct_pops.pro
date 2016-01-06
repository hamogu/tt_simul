; Replacement only for numerical tests!!
; compile to replace original CHIANTI routine, if TT_simul construct of ioneq fails, disables level pop corrections!!! 

FUNCTION correct_pops, pp, t, xne, ionrec, cc, crate=crate, rrate=rrate, $
                       correction=correction, frac_low=frac_low, $
                       frac_high=frac_high
		      
common tt_correct_pops, warning_displayed

if n_elements(warning_displayed) eq 0 then warning_displayed=0

if warning_displayed ne 1 then begin
  print, "USING REPLACEMENT FOR CORRECT POPS: INACCURATRE PHYSICS!!! ONLY FOR TEST PURPOSES!!!"
  warning_displayed=1
endif  

RETURN, pp
end