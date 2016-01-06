; temp=[.1,.2,.3,.4,.6,.8,1.,1.2,1.5,1.8,2.1,2.5,3.,4.,5.,6.,7.5,9.,11.]
; makechiantitable,temp,'test.fits'


pro makechiantitable,intpar1,fitsfile,debug=debug
common elements, abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

; intpar2=[2,3,4,5,6,7,8]
naddparam=8
; fitsfile='name.fits'
addnames=['C','N','O','Ne','Mg','Si','S','Fe']



common tt_units

sngl_ion={master:'/data/hspc62/steh305/idl/tt_results/twhya/master2.ions',C:['C_6','C_5','C_4','C_4d','C_5d','C_3','C_2'],N:['N_6','N_7','N_4','N_3','N_2'],O:['O_7','O_8','O_6','O_6d','O_7d','O_5','O_4','O_3'],Neon:['Ne_9','Ne_10','Ne_8','Ne_8d','Ne_9d','Ne_7','Ne_6','Ne_5'],Mg:['Mg_11','Mg_12','Mg_10','Mg_10d','Mg_9','Mg_8','Mg_7','Mg_6'],Si:['Si_13','Si_14','Si_12','Si_12d','Si_11','Si_10','Si_9','Si_8','Si_7'],S:['S_10','S_11','S_12','S_13','S_14','S_14d','S_9','S_8','S_7','S_6','S_5'],Fe:['Fe_17','Fe_18','Fe_19','Fe_20','Fe_21','Fe_22','Fe_23','Fe_24','Fe_17d','Fe_18d','Fe_19d','Fe_20d','Fe_21d','Fe_22d','Fe_23d','Fe_24d','Fe_16','Fe_15','Fe_14','Fe_13','Fe_12','Fe_11','Fe_16d','Fe_15d','Fe_14d','Fe_13d','Fe_12d','Fe_11d']}
naddparm=8


base_abund='/usr/local/hssoft/chianti/dbase/abundance/version_3/grevesse_sauval98.abund'


;--- generate wavelength and energy scales ---
min_lambda=1.
max_lambda=500.
ang2kev=12.39854
kev_scale=findgen(12000)*0.00025+0.000125
lambda_scale=ang2kev/kev_scale
lambda_scale=lambda_scale[where(lambda_scale ge min_lambda and lambda_scale le max_lambda)]
n_energies=11008

binboundaries=mid2bound(lambda_scale)
upper=binboundaries[0:n_elements(lambda_scale)-1]
lower=binboundaries[1:*]
bin_width=upper-lower
e_max=ang2kev/lower	;lower wavelength->higher energy!
n_energies=n_elements(e_max)
e_min=[0,e_max[0:n_energies-2]]

;--- Setup parameters

nintparm=1
intpoints=n_elements(intpar1)


;--- Fits header ---
fxhmake,header,/initi,/extend,/date
fxaddpar, header, 'MODLUNIT','photons/cm^2/s','Unit of spectrum'
fxaddpar,header,'MODLNAME','VCHIANTI'
fxaddpar,header,'ORIGIN','CHIANTI','by Moritz Guenther, Hamburger Sternwarte'
fxaddpar,header,'REDSHIFT','F','redshift is not a fit parameter'
fxaddpar,header,'ADDMODEL','T','additive model'
fxaddpar,header,'HDUCLASS', 'OGIP    ','format conforms to OGIP standard'
fxaddpar,header,'HDUDOC  ', 'OGIP/92-009','document defining format'
fxaddpar,header,'HDUCLAS1', 'XSPEC TABLE MODEL','model spectra for XSPEC'
fxaddpar,header,'HDUVERS1', '1.0.0   ','version of format'

;--- Parameters section ---
fxaddpar,header1,'EXTNAME','PARAMETERS'
fxaddpar,header1,'NINTPARM',1
fxaddpar,header1,'NADDPARM',naddparm
fxaddpar,header1,'HDUCLASS', 'OGIP    ','format conforms to OGIP standard'
fxaddpar,header1,'HDUCLAS1', 'XSPEC TABLE MODEL','model spectra for XSPEC'
fxaddpar,header1,'HDUCLAS2', 'PARAMETERS','extension containing parameter info'
fxaddpar,header1,'HDUVERS1', '1.0.0   ','version of format'

parameters=replicate({NAME:'Temp          ',METHOD:long(0),INITIAL:1.,DELTA:.1,MINIMUM:min(intpar1),BOTTOM:min(intpar1),TOP:max(intpar1),MAXIMUM:max(intpar1),NUMBVALS:long(n_elements(intpar1)),VALUE:intpar1},naddparm+nintparm)

; elements
for i=nintparm,naddparm+nintparm-1 do begin
  parameters[i].name=addnames[i-nintparm]
  parameters[i].initial=1.
  parameters[i].delta=1.
  parameters[i].minimum=0.
  parameters[i].bottom=.1
  parameters[i].top=3.
  parameters[i].maximum=10.
  parameters[i].numbvals=long(1)
endfor

;--- Energy section ---
fxaddpar,header2,'EXTNAME','ENERGIES'  
fxaddpar,header1,'HDUCLASS', 'OGIP    ','format conforms to OGIP standard'
fxaddpar,header1,'HDUCLAS1', 'XSPEC TABLE MODEL','model spectra for XSPEC'
fxaddpar,header1,'HDUCLAS2', 'ENERGIES','extension containing energy bin info
fxaddpar,header1,'HDUVERS1', '1.0.0   ','version of format'
energies=replicate({ENERG_LO:0.,ENERG_HI:0.},n_energies)


;--- Spectra section ---
fxaddpar,header3,'EXTNAME','SPECTRA'
fxaddpar,header1,'HDUCLASS', 'OGIP    ','format conforms to OGIP standard'
fxaddpar,header1,'HDUCLAS1', 'XSPEC TABLE MODEL','model spectra for XSPEC'
fxaddpar,header1,'HDUCLAS2', 'MODEL SPECTRA','extension containing model spectra'
fxaddpar,header1,'HDUVERS1', '1.0.0   ','version of format'

basespectra={PARAMVAL:[0.,0.],INTPSPEC:fltarr(n_energies)}
for i=1,naddparm do basespectra=add_tag(basespectra,fltarr(n_energies),'ADDSP00'+strtrim(string(i),2))

spectra=replicate(basespectra,intpoints)


; The spectra are given in photons/s/Ang/sr/cm^2, but XSPEC expects photons/s/cm^2
debug=make_array(naddparam+1,n_energies)
element=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si', $
         'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co', $
         'Ni','Cu','Zn']
read_abund,base_abund,abund,abund_ref
; --- Write new abund file ---
openw,af,'temp_master.abund',/get_lun
for i=0,29 do if abund[i] gt 0. and where(strupcase(addnames) EQ strupcase(element[i])) eq -1 then printf,af,i+1,alog10(abund[i])+12.,element[i], format='(i4,f8.2,a6)'
printf,af,'-1'
printf,af,'%file:' ,'temp_master.abund'
printf,af,'%intermedient file for generation of single ion continua; Moritz Guenther'
printf,af,'% file can be deleted'
printf,af,'-1'
free_lun,af


factor=1e+14*4.*!pi*bin_width

for ipar1=0,n_elements(intpar1)-1 do begin
  spectra[ipar1].paramval=intpar1[ipar1]
  tt_synthetic, min_lambda,max_lambda,logem=0,logt_i=alog10(intpar1[ipar1]/8.617e-8), /phot,rphot=1,radtem=6e3,output=lines, master=sngl_ion.master,dens=1e10,ioneq_n='/usr/local/hssoft/chianti/dbase/ioneq/mazzotta_etal_ext.ioneq'
  make_chianti_spec, lines,reverse(lambda_scale),spectrum,/cont,/phot,abund='temp_master.abund'
  spectra[ipar1].intpspec=factor*reverse(spectrum.spectrum)
  debug[0,*]=reverse(spectrum.continuum)
  for k=0,naddparm-1 do begin
    tt_synthetic, min_lambda,max_lambda,logem=0,logt_i=alog10(intpar1[ipar1]/8.617e-8), /phot,rphot=1,radtem=6e3,output=lines,sngl=sngl_ion.(k+1),dens=1e10,ioneq_n='/usr/local/hssoft/chianti/dbase/ioneq/mazzotta_etal_ext.ioneq'
    read_abund,base_abund,abund,abund_ref
    iz=where(strupcase(addnames[k]) EQ strupcase(element))
    ; --- Write new abund file ---
    openw,af,'temp_elem.abund',/get_lun
    printf,af,iz+1,alog10(abund[iz])+12.,element[iz], format='(i4,f8.2,a6)'
    printf,af,'-1'
    printf,af,'%file:' ,'temp_elem.abund'
    printf,af,'%intermediate file for generation of single ion continua; Moritz Guenther'
    printf,af,'% file can be deleted'
    printf,af,'-1'
    free_lun,af

    make_chianti_spec, lines,reverse(lambda_scale),spectrum,abund='temp_elem.abund',/phot,/cont
    spectra[ipar1].(k+2)=factor*reverse(spectrum.spectrum)*1e-12
    debug[k+1,*]=reverse(spectrum.continuum)*1e-12
  endfor
endfor

;file_delete,'temp_elem.abund','temp_master.abund'

for i=0,n_energies-1 do energies[i].energ_lo=e_min[i]
for i=0,n_energies-1 do energies[i].energ_hi=e_max[i]


mwrfits,test,fitsfile,header,/create
mwrfits,parameters,fitsfile,header1
mwrfits,energies,fitsfile,header2
mwrfits,spectra,fitsfile,header3

end