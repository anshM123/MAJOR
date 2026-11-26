#!/usr/bin/env python3
"""Certified spectral-degree cutoff for both channels, conditional on retained inputs.
The proof is in SPECTRAL_DEGREE.md. All numerical inequalities here are exact rational.
"""
from fractions import Fraction as F
from math import factorial,log2
from resid_budget import pow2exp
from base_bounds import logl
from tight_spectral_tail import tail
import json
M=6000;m=M+1;V=F(1030);tau=F(3150);beta=F(19,5);r=F(91,100);c=F(9,10);H=F(1,998**2);sig=F('.901');s2=sig**2
assert sum(beta**(2*j)/factorial(2*j) for j in range(9))>(1+r)/(1-r)
assert (1200**2+F(1,4))<tau*tau
assert F(1,4)+F(49,4)<V*V
# majorant of the angular coefficient operator on the raw f-series, and ideal series
x=c/r
Craw=(1+20*H)*4**6*factorial(6)/r**3/(1-x)**7+H*5**8*factorial(8)/r**4/(1-x)**9
Cideal=H/(1-r)**4
C=Craw+Cideal
assert C<10**25
# Absolute spectral moment (Gaussian and log parts of perturbed raw weight;
# ideal G mass <= 6e7 and imaginary-frequency atom total mass <= 1).
t=2*m/V
log_exponent=(t+8)**2/2-t*(V-1)
assert log_exponent<0
Mass=2*F('2.27')*1220001*pow2exp(s2*t*t/2)+30+60000000+1
near=pow2exp(beta*tau)*(V/tau)**(2*m)*Mass
# The unrestricted positive majorant for real frequencies beyond 1200.
a=F(1200)
q=beta*(a+1)-(a-1029)**2/(2*s2);sl=(a-1029)/s2-beta
q2=beta*(a+1)+8*a-a*a/2;sl2=a-beta-8
assert min(sl,sl2)>0
far=2*1220001*pow2exp(q)/sl+10*pow2exp(q2)/sl2
E=C*(near+far)
assert E<F(1,2**1910)
# For the actual (unperturbed) spectrum, ordinary-power truncation is still charged.
assert tail(c,17199)<F(1,2**140)
print('PASS_SPECTRAL_DEGREE',M,'raw operator bound',float(Craw),'log2 error upper',log2(E.numerator)-log2(E.denominator),'ordinary K',17199,'ordinary tail log2',log2(tail(c,17199).numerator)-log2(tail(c,17199).denominator))
with open('spectral_cut_results.json','w')as f:json.dump({'M':M,'ordinary_K':17199,'Craw':str(Craw),'Cideal':str(Cideal),'spectral_error_upper_log2_rounded_up':__import__('math').ceil(log2(E.numerator)-log2(E.denominator)),'status':'PASS_CONDITIONAL_SPECTRAL_DEGREE_ONLY'},f,indent=2)
