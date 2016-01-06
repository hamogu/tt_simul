;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	PREPARE_ANALYSIS
;
; PURPOSE:
;	This procedure calculates the basic diagnostcs of a plasma from a simulated hydrodynamical structure.
;
; CATEGORY:
;	analysis
;
; CALLING SEQUENCE:
;	prepare_analysis, hydrodyn, delta_x,volem,elecdens,o7i,o7f,o7r,o7f2i,ne9i,ne9f,ne9r,ne9f2i,quiet=quiet,radtemp=radtemp,rphot=rphot
;
; INPUTS:
;	hydrodyn:	hydrodynmic shock structure as output by tt_simul
;
; OPTIONAL INPUTS:
;	radtemp:	temperatue of a background black-body radiation field [Kelvin]
;	rphot:		dilution factor, describing the dilutes of the background radiation field [1 is stellar surface]
;
; OUTPUTSS:
;	delta_x:	a vector with the stepsize between different point in hydrodyn
;	volem:		a vector containing the volume emission meassure of every step [cm^{-3}]
;	elecdens:	a vector containing the electron density in every step [cm^{-3}]
;	o7i,o7f,o7r,ne9i,ne9f,ne9r:	vectors containig the emissivity in every step [erg/s/cm^3] in the lines marked
;	o7f2i,ne9f2i:	vectors containig the f/i ratio for every step
;
; KEYWORDS:
;	quiet:		suppress printout
;
; TEST STATUS:
;	checked
;
; MODIFICATION HISTORY:
; 	Written by:	Moritz Guenther 22.08.2005
;	hydrodynamics souced out to prepare_analyses_hyd:	3.11.2005  Moritz Guenther
;-

pro prepare_analysis, hydrodyn, delta_x,volem,elecdens,o7i,o7f,o7r,o7f2i,ne9i,ne9f,ne9r,ne9f2i,quiet=quiet,radtemp=radtemp,rphot=rphot
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
prepare_analysis_hyd, hydrodyn, delta_x,volem,elecdens
; O VII
o7lines=calcmygoft(8,7,[[22.,22.2],[21.7,21.9],[21.5,21.7]],reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot)
o7i=o7lines[*,1]
o7f=o7lines[*,0]
o7r=o7lines[*,2]
o7f2i=o7f/o7i
;effective o7f2i, avareged over the whole area
;the first line sometimes contains NaNs, so forget about it
if not keyword_set(quiet) then print, 'f/i that should be observed (O VII)',total(o7f*volem)/total(o7i*volem)
;Ne IX
ne9lines=calcmygoft(10,9,[[13.69,13.7],[13.55,13.57],[13.44,13.45]],reform(elecdens),reform(hydrodyn[5,*]),/quiet,radtemp=radtemp,rphot=rphot)
ne9i= ne9lines[*,1]
ne9f= ne9lines[*,0]
ne9r= ne9lines[*,2]
ne9f2i=ne9f/ne9i
;effective o7f2i, avareged over the whole area
;the first line sometimes contains NaNs, so forget about it
if not keyword_set(quiet) then print, 'f/i that should be observed (Ne IX)',total(ne9f*volem)/total(ne9i*volem)
end