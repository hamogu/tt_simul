pro setplotps,filename
set_plot,'ps'
device,filename=filename+'.ps',encapsulated=0
end