pro setploteps,filename
set_plot,'ps'
device,filename=filename+'.eps',/encapsulated
end