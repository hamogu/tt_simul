pmultiold=!p.multi
!p.multi=[0,1,1]
!p.charsize=1.
xmargin=[8,0]
ymargin=[2,2]
xychar=1.3
xl=0.1
xr=0.1
yu=0.1
yo=0.1
dy=(1.-yu-yo)/3.
dx=(1.-xr-xl)/3.

plot, hyd225[0,*]/hyd225[0,0],ion225[*,9,6],yr=[0,1],xr=[0,1],pos=[xl+0*dx,yu+2*dy,xl+1*dx,yu+3*dy],charsize=1e-4
axis, yaxis=0,ytit='rel abundance'
xyouts,xl+.2*dx,yu+3.1*dy,textoidl('v_0=225 km/s'),/normal
xyouts,xl+1.2*dx,yu+3.1*dy,textoidl('v_0=350 km/s'),/normal
xyouts,xl+2.2*dx,yu+3.1*dy,textoidl('v_0=500 km/s'),/normal
oplot, hyd225[0,*]/hyd225[0,0],ion225[*,9,6]+ion225[*,9,7]
oplot, hyd225[0,*]/hyd225[0,0],ion225[*,9,6]+ion225[*,9,7]+ion225[*,9,8]
xyouts, 0.3,0.05, 'Ne VII',chars=xychar
xyouts,0.3,.4,'Ne VIII',chars=xychar
xyouts, .3,.7,'Ne IX',chars=xychar
plot, hyd350[0,*]/hyd350[0,0],ion350[*,9,8],tit=textoidl('v_0=350 km/s'),xr=[0,1],pos=[xl+1*dx,yu+2*dy,xl+2*dx,yu+3*dy],/noerase,chars=1e-4
xyouts,.3,.4,'Ne IX',chars=xychar
plot, hyd500[0,*]/hyd500[0,0],ion500[*,9,8],tit=textoidl('v_0=500 km/s'),xr=[0,1],pos=[xl+2*dx,yu+2*dy,xl+3*dx,yu+3*dy],/noerase,chars=1e-4
oplot, hyd500[0,*]/hyd500[0,0],ion500[*,9,8]+ion500[*,9,9]
xyouts,.2,.3,'Ne IX',chars=xychar
xyouts,.2,.78,'Ne X',chars=xychar
plot, hyd225[0,*]/hyd225[0,0],ne9i225/max(ne9i225),xr=[0,1],pos=[xl+0*dx,yu+1*dy,xl+1*dx,yu+2*dy],/noerase,charsize=1e-4
;axis,yaxis=0,ytit=textoidl('rel vol emissivity')
xyouts, .25*xl,yu+dy,textoidl('rel vol emissivity'),orient=90,/normal
xyouts, .55*xl,yu+dy,textoidl(' [arbitrary units]'),orient=90,/normal
legend,['f','i'],/center,linestyle=[0,1],pos=[.4,.9]
oplot, hyd225[0,*]/hyd225[0,0],ne9f225/max(ne9i225),line=1
plot, hyd350[0,*]/hyd350[0,0],ne9i350/max(ne9i350),xr=[0,1],pos=[xl+1*dx,yu+1*dy,xl+2*dx,yu+2*dy],/noerase,chars=1e-4
oplot, hyd350[0,*]/hyd350[0,0],ne9f350/max(ne9i350),line=1
plot, hyd500[0,*]/hyd500[0,0],ne9i500/max(ne9i500),xr=[0,1],pos=[xl+2*dx,yu+1*dy,xl+3*dx,yu+2*dy],/noerase,chars=1e-4
oplot, hyd500[0,*]/hyd500[0,0],ne9f500/max(ne9i500),line=1
yyplot, hyd225[0,*]/hyd225[0,0],hyd225[5,*],hyd225[2,*],/yllog,ylrange=[1e5,1e7],yrrange=[1e12,2e13],xr=[0,1],yltit='Temperature [K]',pos=[xl+0*dx,yu+0*dy,xl+1*dx,yu+1*dy],/noerase,rchar=1e-4
legend,['temp','density'],/top,linestyle=[0,2]
yyplot, hyd350[0,*]/hyd350[0,0],hyd350[5,*],hyd350[2,*],/yllog,ylrange=[1e5,1e7],yrrange=[1e12,2e13],xr=[0,1],pos=[xl+1*dx,yu+0*dy,xl+2*dx,yu+1*dy],/noerase,chars=1e-4
axis,xaxis=0,xtickv=[.2,.4,.6,.8,1.0],xticks=4
yyplot, hyd500[0,*]/hyd500[0,0],hyd500[5,*],hyd500[2,*],/yllog,ylrange=[1e5,1e7],yrrange=[1e12,2e13],xr=[0,1],yrtit=textoidl('density [cm^{-3}]'),pos=[xl+2*dx,yu+0*dy,xl+3*dx,yu+1*dy],/noerase,lchars=1e-4
axis,xaxis=0,xtickv=[.2,.4,.6,.8,1.0],xticks=4
!p.multi=pmultiold                                            
delvar, pmultiold