openw,1,'3d.data'
printf,1,'# Visualising maxdepth in 3D - data for use with gnuplot'
for i=0, 16 do begin & for j=0,12 do printf,1, (sort2d(rad6e3.maxdepth,v0,n0))[i,j]/km,(sort2d(v0,v0,n0))[i,j],(sort2d(n0,v0,n0))[i,j] & printf,1,'' & endfor  
close,1
