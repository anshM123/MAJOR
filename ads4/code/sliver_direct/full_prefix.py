#!/usr/bin/env python3
"""Uniform physical prefix truncation over all remaining L below 2.3.
No claim of sign. Uses retained spectral envelopes and angular identity.
"""
from fractions import Fraction as F
from base_bounds import sqrtlo,logl,hi
from resid_budget import gaussian_mass
from math import ceil
import json
R=F(1999,2000);DEL=F(1,4000);V=1100;AMIN=998;SMAX=9
# radius |z|=R and Cauchy dressing at radius DEL
assert R+DEL==F(3999,4000)
assert 1+(2/DEL**2+5/DEL+8)/AMIN**2<34
# H uses d^3/(1+qh)^4<=1; numerator dressing qh <=20/998^2.
# use generous 40 times the undressed coefficient bound for dressing.
# 2/Lambda integral |f_theta| <60000001, tau Lambda^2 integral |psi3|<60000000.
assert 80*4000**4*120000001 < 2**83
assert F(22,7)*V/2<2494*logl(F(2))
# combined coefficient envelope 2^(2577) on radius R
POW=2700
assert 1/(1-F('.9')/R)<11
elo=F(11,2**(POW-2577))
ehi=880*4000**4*gaussian_mass(V,F(11,7))
assert elo+ehi<F(1,2**110)
# Check L range: exp(-2.3) >0.1 and exp(-.06)<.95
from check_lowfreq import exp_u
assert 1/exp_u(F('2.3'))>F('.1')
assert exp_u(-F('.06'))<F('.95')
def need_k(c):
    t=c/R;p,q=t.numerator,t.denominator
    def ok(k):return (p**(k+1)<<POW)<=q**(k+1)
    lo,up=-1,1
    while not ok(up):up*=2
    while up-lo>1:
        m=(lo+up)//2
        if ok(m):up=m
        else:lo=m
    return up

def cap(c):return hi((1-sqrtlo(c))**2/AMIN**2)
def gap(l,u):
    K=need_k(u);h=cap(l)
    # denominator includes through index K-1, matching U_K
    return K,1-h*(K-1)*(K-1+SMAX),h
c=F('.05');end=F('.90');rows=[]
while c<end:
    l=0;u=int((end-c)*10000)+1
    # monotone in c_hi at fixed c_lo
    while u-l>1:
        j=(l+u)//2;v=c+F(j,10000)
        k,g,h=gap(c,v)
        if g>=F(1,20):l=j
        else:u=j
    assert l>0
    v=c+F(l,10000);k,g,h=gap(c,v)
    row=dict(c_lo=str(c),c_hi=str(v),K=k,h_cap=str(h),gap=str(g),gap_approx=float(g));rows.append(row);print(json.dumps(row));c=v
assert c==end
print('PASS_PREFIX_ONLY:',len(rows),'cells; Kmax',max(x['K'] for x in rows),'min_gap',min(x['gap_approx'] for x in rows))
json.dump(dict(status='PASS_PREFIX_ONLY',c_lo='.05',c_hi='.9',L_includes=['.06','2.3'],R=str(R),pow=POW,error='2^-110',rows=rows),open('full_prefix_results.json','w'),indent=2)
