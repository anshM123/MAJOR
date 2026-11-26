#!/usr/bin/env python3
"""Uniform physical finite-prefix reduction for the remaining small-L band.
No assumption on a finite-mass h-series or connection-series localisation.
Uses the retained angular identity, real spectral residual/cutoff and Gaussian envelopes.
This certifies truncation and absence of poles, NOT the sign of the prefix.
"""
from fractions import Fraction as F
from base_bounds import sqrtlo,logl,hi
from resid_budget import gaussian_mass
import json
R=F(99,100);DEL=F(1,200);V=1100;AMIN=998;SMAX=9
POW=2670
# F_{p,nu} and derivatives on |z|=R, using the Mehler estimate on |z|<=R+DEL.
assert R+DEL==F(199,200)
assert 1+(2/DEL**2+5/DEL+8)/AMIN**2<2
# (2/Lambda) int |f_theta| < 60000001, tau*Lambda^2 int |psi3|<60000000.
# 2*Mehler prefactor (pi/2) <4, (1-R-DEL)^(-4)=200^4.
assert 4*200**4*120000001 < 2**60
# e^(pi*1100/2) < 2^2494 (pi<22/7 and lower bound for log 2).
assert F(22,7)*V/2<2494*logl(F(2))
assert 1/(1-F('.56')/R)<3
# Physical U_k <= c^k, including the zero beyond termination.
# Thus low-frequency error < 3*2^(2554-2670), high-frequency error below.
elo=F(3,2**(POW-2554))
# For the tail above 1100 use whole kernel/prefix bound by geometric series.
# 4*200^4, 2/Lambda<1, geometric factor<3. Only k2 contributes there.
ehi=12*200**4*gaussian_mass(V,F(11,7))
assert elo+ehi<F(1,2**110)

def need_k(c):
    # smallest K with (c/R)^(K+1) <= 2^-POW; exact integer powers
    t=c/R; p,q=t.numerator,t.denominator
    def ok(k):return (p**(k+1)<<POW)<=q**(k+1)
    lo,up=0,1
    while not ok(up):up*=2
    while up-lo>1:
        m=(lo+up)//2
        if ok(m):up=m
        else:lo=m
    return up

def cap(c):return hi((1-sqrtlo(c))**2/AMIN**2)
def pole_gap(lo,hi_):
    K=need_k(hi_);h=cap(lo);return K,1-h*(K-1)*(K-1+SMAX),h
# Simple coarse partition; maximize c_hi subject to a certified denominator gap.
c=F('.05');end=F('.56');rows=[]
while c<end:
    u=c+F('.01');best=None
    while u<=end:
        K,g,h=pole_gap(c,u)
        if g<F(1,10):break
        best=(u,K,g,h);u+=F('.01')
    assert best is not None
    u,K,g,h=best
    row=dict(c_lo=str(c),c_hi=str(u),K=K,h_cap=str(h),denominator_gap=str(g),denominator_gap_approx=float(g))
    rows.append(row);print(json.dumps(row));c=u
assert c==end
out=dict(c_lo='0.05',c_hi='0.56',nu_cut=V,coefficient_circle=str(R),prefix_error_lt='2^-110',low_error_upper=str(elo),high_error_log2_upper=-(ehi.denominator.bit_length()-ehi.numerator.bit_length()-1),rows=rows,status='TRUNCATION_AND_POLE_COVER_ONLY')
open('low_prefix_results.json','w').write(json.dumps(out,indent=2)+'\n')
print('PASS truncation and pole exclusion:',len(rows),'cells; maximum K',max(x['K'] for x in rows),'minimum gap',min(x['denominator_gap_approx'] for x in rows))
