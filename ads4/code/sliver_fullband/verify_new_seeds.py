#!/usr/bin/env python3
"""Exact input, fixed-point range and coherent Gaussian seed-error checks."""
from fractions import Fraction as F
from math import ceil
import json,sys,glob
from resid_budget import pow2exp
sys.set_int_max_str_digits(0)
w=json.load(open('weights.json'));cs=w['centres_sigma'];xs=w['x_re'];ys=w['y_re'];assert len(cs)==len(xs)==3564
assert all(0<=F(c)<=1029 and F('.449')<F(s)<F('.901')for c,s in cs)
assert all(abs(F(y))<1 for y in ys)
W=sum(abs(F(x))for x in xs);assert W<1220000
sig=F('.901');B=F(11,7);K=6000;D=3000
# Two boundary sequences and the main normalized moment recurrence, up to 12000.
# The boundary mean is a convex combination of the center and +/-nu_c, and |nu_c|<1099.
assert 997002<999**2
boundary_a=F(1099,1200)+2*K*sig**2/1200**2
main_a=F(1029,1200)+2*K*sig**2/1200**2
b=sig**2/1200
assert boundary_a<F('.93') and main_a<F('.87') and b<F('.001')
# |M0|<=2.27+1e-780, main true sequence magnitude <=3, boundary true sequence <=1.
assert 2*F(22,7)*sig**2<F('2.27')**2
assert main_a*3+2*b<3
# A roundoff in a boundary mean costs <= 2/S times |approx moment|<=2.
# Initial rounded boundary seeds have <=2/S error. Both products truncate once.
assert 100*boundary_a+4+2<100
assert 100*main_a+200*b+3<100
# Final weighted summation allowance per shard.
assert W*100+3564<10**10
seeds=sorted(glob.glob('fixed_6000_3000_*.seeds'));assert len(seeds)==12
n=0
for name in seeds:
    for line in open(name):
        l,u=map(int,line.split());assert 0<=l<=u and u-l<=10**20;n+=1
assert n==3564
# Every chosen midpoint for M0 is within 1e-780 of the exact seed. Its homogeneous
# error is an ordinary Gaussian moment sequence (normalization mass delta).
DM=W*F('1e-780');realW=2*DM*pow2exp(B*1029+B*B*sig**2/2)
# Inherited real-frequency raw action kernel bound, independently proved in
# raw_seed_budget_full.py; use a greatly inflated 1e19 and factor 2 for qh.
C=2*10**19
coeff={1:F(-96),9:F(5184,5),49:F(-4704,5)}
realErr=C*realW*sum(abs(a)/p for p,a in coeff.items())/998
# Keeping the original rational-pole seed fixed also creates a pole atom at each pole.
# Its mass is at most DM/a (integrate positive Gaussian against 1/(4nu^2+a)).
rr=json.load(open('pole_bound_full_results.json'))
for b in (0,4):
    poleErr=sum(abs(coeff[x['a']])*DM/x['a']*F(x['action_bound'])/998 for x in rr if x['b']==b)
    assert realErr+poleErr<F('1e-45')
    print('PASS Gaussian M0 coherent error b',b,'<1e-45; real bound',float(realErr),'pole bound log10',__import__('math').log10(poleErr.numerator)-__import__('math').log10(poleErr.denominator))
print('PASS all 3564 seeds and fixed-point constants; 12 shards')
