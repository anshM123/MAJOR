#!/usr/bin/env python3
"""Rigorous checks for Lemma G1 (physical twist cap).
(1) Symbolic: (2x+s)^2 - 4x(x+s) == s^2  (so f(x)=sqrt(x(x+s)) has f'(x) >= 1 for x>0, s>=0).
(2) Exact-rational stress test of the conclusion on integer/rational states (not needed for the proof,
    a guard against algebra slips): for n in a grid >= 1001-b, a in a rational grid >= 998:
       h <= c/((1001-b)(1002+b))   and   998^2 h <= (1-sqrt(c))^2,
    the latter checked exactly as: sqrt(c) <= 1-998 sqrt(h)  <=>  sqrt(n(n+s)) + 998 <= sqrt(N(N+s))
    <=> (sqrt(n(n+s))+998)^2 <= N(N+s)  <=> 2*998*sqrt(n(n+s)) <= N(N+s)-n(n+s)-998^2 (rhs>0 checked)
    <=> 4*998^2*n(n+s) <= (N(N+s)-n(n+s)-998^2)^2   (all exact Fractions).
(3) The driver rounding: s_lo = isqrt(floor(c_lo 10^40))/10^20 satisfies s_lo^2 <= c_lo, hence
    (1-s_lo)^2 >= (1-sqrt(c_lo))^2; cap rounded up. Checked for every cell endpoint used (imported from driver).
"""
import sympy as sp
from fractions import Fraction as F
import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from gap_triangle import g1_cap
x,s=sp.symbols('x s',positive=True)
assert sp.expand((2*x+s)**2-4*x*(x+s)-s**2)==0
cnt=0
for b in (0,4):
 S=2*b+1
 for n in [1001-b,1002-b,1003,1500,2000,5000,10**4,10**5,10**6,10**8]:
  for a in [F(998),F(9981,10),F(999),F(1000),F(1234567,1000),F(5000),F(10**4),F(10**6),F(10**9)]:
   N=n+a;h=1/(N*(N+S));c=n*(n+S)*h
   assert h<=c/((1001-b)*(1002+b))
   rhs=N*(N+S)-n*(n+S)-998**2;assert rhs>0 and 4*998**2*n*(n+S)<=rhs*rhs
   cnt+=1
for l in ['.28','.29','.3','.35','.4','.45','.5','.51','.55','.56']:
 for b in (0,4):
  hc,c1,c2=g1_cap(F(l),F(l)+F(1,100),b)
  # exact: c2 >= (1-sqrt(l))^2/998^2  <=>  1 - 998 sqrt(c2) <= sqrt(l)  (check by squaring when lhs>=0)
  t=F(l);u=c2*998**2  # need sqrt(u) >= 1-sqrt(t), i.e. sqrt(t) >= 1-sqrt(u)
  # verify with rational bracket: find r with r^2<=u exactly and (1-r)^2<=t, r<=1
  from math import isqrt
  r=F(isqrt((u.numerator*10**60)//u.denominator),10**30);assert r*r<=u and r<1 and (1-r)**2<=t
  cnt+=1
print('G1 checks PASS:',cnt,'exact instances; symbolic derivative identity OK')
