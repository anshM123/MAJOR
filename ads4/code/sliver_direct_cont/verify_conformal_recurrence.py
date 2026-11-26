#!/usr/bin/env python3
"""Symbolic recurrence derivation and exact finite-spin sanity checks.
Not a positivity certificate for the remaining interval.
"""
from fractions import Fraction as F
from math import comb
import sympy as S
w,j,c,h,s=S.symbols('w j c h s');z=4*w/(1+w)**2
p=(1-w*w)**2;fp=j*p+w*S.diff(p,w);fpp=(j-1)*fp+w*S.diff(fp,w)
cs=[-c+h*j*(j+s),4-2*h*j*(s-1),2*c-2*h*j*j-8*h*j-4*h*s-8*h-8,4+2*h*(j+4)*(s-1),-c+h*(j+4)*(j+4-s)]
assert S.cancel((z-c)*p+h*(fpp+(1+2*w/(1-w*w)+s*(1-w)/(1+w))*fp)-sum(x*w**k for k,x in enumerate(cs)))==0
checks=0
for s0 in [1,9]:
 for n in range(1,31):
  for a in [2,998]:
   h0=F(1,(n+a)*(n+a+s0));c0=n*(n+s0)*h0
   U=[F(1)]
   for k in range(n):U.append(U[-1]*(c0-h0*k*(k+s0))/(1-h0*k*(k+s0)))
   V=[F(1)]+[sum(F(j0*comb(2*k,k-j0),k*4**k)*U[k] for k in range(j0,n+1)) for j0 in range(1,n+1)]+[F(0)]*4
   assert all(x>=0 for x in V)
   for j0 in range(n+1):
    co=[-c0+h0*j0*(j0+s0),4-2*h0*j0*(s0-1),2*c0-2*h0*j0*j0-8*h0*j0-4*h0*s0-8*h0-8,4+2*h0*(j0+4)*(s0-1),-c0+h0*(j0+4)*(j0+4-s0)]
    assert sum(q*V[j0+k] for k,q in enumerate(co))==0;checks+=1
print('PASS: symbolic identity and',checks,'exact finite-spin recurrence checks')
