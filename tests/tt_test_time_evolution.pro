;+
; PROJECT:
;	TT_SIMUL
; NAME:
;	TT_TEST_TIME_EVOLUTION
;
; PURPOSE:
;	TT_test_ioneq and TT_test_ioneq_rate check if the ionisation stage calculation relaxes to its equlilibrium value for long times.
;	If it does the rates used are correct. This procedure checks that the time scales are realistic.
;
; CATEGORY:
;	VERIFICATION
;
; CALLING SEQUENCE:
;	tt_test_time_evloution,temp, startH1,timescale,nH1
;
;  INPUTS:
;	temp:	Teperature in K
;	startH1:	number between 0. and 1., is ioneq[0,0,0] for time=0
;			ioneq[0,0,1]=1-startH1;      
; OUTPUTS:
;	timescale:	typical timescale
;	nH1:		array with time evolution for 4*timescale in 0.1 steps
;
; OPTIONAL OUTPUTS:
;	n1inf:		the equilibrium value for HI: n(HI,t=inf)
;                        
; EXAMPLE:
;	tt_test_time_evloution,10.^4.2, 0.,timescale,nh1
;
;  COMNON BLOCKS:
;	tt_steps, elements, tt_mean_atomic_weight
;
; CALLING:
;	tt_dielec_recomb,tt_rad_recomb,tt_arnaudraymond,tt_get_mean_atomic_weight, tt_microscopic
;
; TEST STATUS:  
;	tested and used
;
; MODIFICATION HISTORY:
;	Written by:	Moritz Günther 28.4.2005
;	Moritz Günther 21.6.2005:	added n2inf output
;-
pro tt_test_time_evolution,temp, startH1,timescale,nH1,n1inf=n1inf
common tt_steps, delta_x, A, delta_t,x
common elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref
common tt_mean_atomic_weight, mean_atomic_weight

x=0.
A=1.
delta_x=1.

n=1e11

abund=make_array(50,val=0.)
abund[0]=1.			;pure H
abund_ref='Test case for time evolution'

ioneq=make_array(1,30,31)
ioneq[0,0,0:1]=[startH1,1.-startH1]
ioneq_logt=alog10(temp)
ioneq_ref='Test case for time evolution'

mean_atomic_weight=tt_get_mean_atomic_weight()
x_e=tt_get_xe(Temp)
out=obj_new('tt_output','../tt_results/test_rate','Test time evolution 1.0')

; --- analytic for H ---
recomb=make_array(31,31,val=0.)
colion=make_array(31,31,val=0.)
;there is no dielectronoc recombination for H
tt_rad_recomb,Temp,recomb
tt_arnaudraymond, Temp,colion
R12=colion[0,0]
R21=recomb[0,1]
timescale=1./(R12*n)
n1inf=R12/(R12+R21)

delta_t=timescale/10.

nH1=make_array(41)

nH1[0]=ioneq[0,0,0]
for i=1,40 do begin
  tt_microscopic, Temp,Temp,n,x_e,d_xe,en_loss,relabundHII,out,/rate_test
  nH1[i]=ioneq[0,0,0]
endfor  
end