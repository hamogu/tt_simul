
pro testmakechiantitable,intpar1,fitsfile
; intpar2=[2,3,4,5,6,7,8]
naddparm=0



common tt_units



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


;--- Parameters section ---
fxaddpar,header1,'EXTNAME','PARAMETERS'
fxaddpar,header1,'NINTPARM',1
fxaddpar,header1,'NADDPARM',naddparm

parameters={NAME:'Temp          ',METHOD:long(0),INITIAL:1.,DELTA:.1,MINIMUM:min(intpar1),BOTTOM:min(intpar1),TOP:max(intpar1),MAXIMUM:max(intpar1),NUMBVALS:long(n_elements(intpar1)),VALUE:intpar1}


;--- Energy section ---
fxaddpar,header2,'EXTNAME','ENERGIES'  
energies=replicate({ENERG_LO:0.,ENERG_HI:0.},n_energies)


;--- Spectra section ---
fxaddpar,header3,'EXTNAME','SPECTRA'

basespectra={PARAMVAL:[0.],INTPSPEC:fltarr(n_energies)}
spectra=replicate(basespectra,intpoints)




; The spectra are given in photons/s/Ang/sr/cm^2, but XSPEC expects photons/s/cm^2
 factor=1e+14*4.*!pi*bin_width

for ipar1=0,n_elements(intpar1)-1 do begin
  spectra[ipar1].paramval=intpar1[ipar1]
  tt_synthetic, min_lambda,max_lambda,logem=1,logt_i=alog10(intpar1[ipar1]/8.617e-8), /phot,rphot=1,radtem=6e3,output=lines, dens=1e10,ioneq_n='/usr/local/hssoft/chianti/dbase/ioneq/mazzotta_etal_ext.ioneq'
  make_chianti_spec, lines,reverse(lambda_scale),spectrum,/cont,/phot,abund=base_abund
  spectra[ipar1].intpspec=factor*reverse(spectrum.spectrum)
endfor


for i=0,n_energies-1 do energies[i].energ_lo=e_min[i]
for i=0,n_energies-1 do energies[i].energ_hi=e_max[i]


mwrfits,test,fitsfile,header,/create
mwrfits,parameters,fitsfile,header1
mwrfits,energies,fitsfile,header2
mwrfits,spectra,fitsfile,header3

end