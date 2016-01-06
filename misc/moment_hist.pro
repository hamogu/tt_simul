function moment_hist, x,values
mean=total(x*values)/total(values)
variance=total((x-mean)^2.*values)/total(values)	; I use 1/N here not 1/(N-1) as I should, but that is OK for large numbers and avoids trouble with small normalisations, if not a real histogramm is given, but e.g. a spectrum in flux units and not photon counts
skewness=total((x-mean)^3.*values)/(variance^(3./2.)*total(values))
return,[mean,variance,skewness]
end

; IDL> help, gaussfit(o6wave,o6flux,a,chisq=chi,nterms=4,measure=o6err)
; <Expression>    DOUBLE    = Array[308]
; IDL> print, a,chi
;    2.5891449e-13       1032.1955      0.35285946   6.1085853e-15
;        1.7879523
; IDL> g= gaussfit(o6wave,o6flux,a,chisq=chi,nterms=4,measure=o6err)
; IDL> oplot, o6wave,g,col=200
; IDL> oplot, o6wave,g,col=8
; IDL> % LICENSE MANAGER: Server reconnect after 3 attempts (3:0 minutes).
; 
; IDL> save, o6wave,o6flux,o6err,file='o61032line.sav'
