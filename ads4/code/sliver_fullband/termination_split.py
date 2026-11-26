#!/usr/bin/env python3
"""Exact nonphysical witness and a termination-preserving finite reduction."""
from fractions import Fraction as F
from math import ceil
from base_bounds import sqrtlo
import json,sys
sys.set_int_max_str_digits(0)
D=100;S=10**D;K=1567;c=F('.35');h=F('2.33e-7');b=4;s=2*b+1;q=b*(b+1)
# This point is in the rejected rectangle, but violates the twist constraint.
r=json.load(open('next_67_250_7_20_4.json'));assert F(r['c_lo'])<=c<=F(r['c_hi']) and 0<=h<=F(r['hcap'])
assert h>(1-sqrtlo(c))**2/998**2 and (c/h).denominator!=1
with open('coeff_long_17200.txt') as f:
 N,DD=map(int,next(f).split());dv=10**(DD-D);assert N>=K+1 and DD>=D
 A=[];G=[]
 for _ in range(K+2):
  a,ar,g,gr=map(int,next(f).split());A.append((a//dv,ceil(F(ar,dv))+1));G.append((g//dv,ceil(F(gr,dv))+1))
v=e=0
for k in range(K,-1,-1):
 a,ar=A[k];an,anr=A[k+1];g,gr=G[k];k1=(k+1)**2;k2=-k*k+(8-4*b)*k+9-2*b*(b-1)+2*q
 t=F(k1*an+k2*a-2*g,2);te=F(k1*anr+abs(k2)*ar+2*gr,2)
 f=a+h*t;fe=ar+h*te
 if k<K:
  p=k*(k+s);ratio=(c-h*p)/(1-h*p);assert abs(ratio)<=1
  f+=ratio*v;fe+=abs(ratio)*e
 v=int(f);e=ceil(fe)+1
norm=(1-c)**3/(1+q*h)**4
lb=F(v-e,S)*norm;ub=F(v+e,S)*norm
assert ub<-10**18
print('PASS exact negative NONPHYSICAL rational extension witness:',float(lb),float(ub))
# On physical states use x=h/c=1/(n(n+s)). For n>=K use a continuous
# x-box 0<=x<=1/(K(K+s)); all factors are nonnegative, gaps >=1-c.
# For n<K restrict to the finite exact n's; a>=998 is equivalent to
# c<=n(n+s)/[(n+998)(n+998+s)]. The sum terminates exactly at n.
rows=json.load(open('full_grid_results.json'))['rows'];out=[];Nall={0:set(),4:set()};count=0
for r in rows:
 l,u=F(r['c_lo']),F(r['c_hi']);K=r['K']
 for b in (0,4):
  s=2*b+1;n0=1001-b;finite=[]
  for n in range(n0,K):
   cmax=F(n*(n+s),(n+998)*(n+998+s));up=min(u,cmax)
   if up<l:continue
   finite.append(n);Nall[b].add(n);count+=1
  # Use an interval of n's because cmax is increasing in n.
  assert not finite or finite==list(range(finite[0],K))
  xcap=min(F(1,K*(K+s)),(1-sqrtlo(l))**2/(998**2*l));gap=1-u*xcap*(K-1)*(K-1+s)
  assert gap>=1-u>0
  out.append(dict(c_lo=str(l),c_hi=str(u),K=K,b=b,continuous_xcap=str(xcap),continuous_pole_gap=str(gap),finite_spin_first=finite[0] if finite else None,finite_spin_last=finite[-1] if finite else None,finite_spin_count=len(finite),sign='NOT_RUN'))
res=dict(status='TERMINATION_PRESERVING_REDUCTION_ONLY',witness=dict(b=4,c=str(c),h=str(h),K=1567,lower=str(lb),upper=str(ub),physical=False),rows=out,continuous_rectangles=len(out),finite_spin_cell_count=count,distinct_n={str(b):len(Nall[b])for b in (0,4)})
json.dump(res,open('termination_split_results.json','w'),indent=2)
print('PASS termination-preserving split:',len(out),'channel rectangles,',count,'finite-spin interval checks; distinct n:',res['distinct_n'])
