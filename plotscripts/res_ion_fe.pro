;suited in this way to produce a good .eps outfile
!p.charsize=2.
plot, hydrodyn[0,*]/km,findmax(reform(ioneq[*,25,*])),xr=[-5,140],xstyle=1,xtit='depth [km]',ytit='most abund. ion'

                                          