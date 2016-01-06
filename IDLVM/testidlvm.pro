pro testidlvm
@~/idl_startup.pro
uv=spectral_series([525],[900.,1800.],d_lambda=.1,abund='../../tt_data/twhya2.abund')
save,uv,filename='testidlvm_data.sav'
end