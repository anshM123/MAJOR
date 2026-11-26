#!/usr/bin/env python3
"""Rigorous lower bound on the compact M1 interval using monotone scalar tables.
Decimal.exp is correctly rounded; adjacent decimals enclose it. Every subsequent
bound (including square roots) uses exact rational/integer arithmetic.
"""
from decimal import Decimal,localcontext
from fractions import Fraction as F
from math import isqrt
import json
S=10**18; G=200; N=256

def expbox(x):
 with localcontext() as c:
  c.prec=60
  z=(Decimal(x.numerator)/Decimal(x.denominator)).exp()
  # x is a terminating decimal on this grid, so the argument is exact.
  return F(z.next_minus()),F(z.next_plus())
def ceil(a):return -((-a.numerator)//a.denominator)
def floor(a):return a.numerator//a.denominator
# Tables at t=j/G, up to beta=12.
sl=[];su=[];al=[];au=[];nl=[];dl=[];du=[]
for j in range(12*G+1):
 t=F(j,G)
 if not j:
  sl.append(0);su.append(0);al.append(S);au.append(S);nl.append(floor(S*F(1,24)));dl.append(S);du.append(S);continue
 el,eu=expbox(t);hl,hu=expbox(t/2)
 # monotone interval operations for sinh/cosh
 sslo=(el-1/el)/2;sshi=(eu-1/eu)/2
 cl=(el+1/eu)/2;cu=(eu+1/el)/2
 hcl=(hl+1/hu)/2
 sl.append(floor(S*sslo));su.append(ceil(S*sshi))
 al.append(floor(S*sslo/t));au.append(ceil(S*sshi/t))
 nl.append(floor(S*(t*t*hcl-2*(cu-1))/t**4))
 dl.append(floor(S*2*(cl-1)/t**2));du.append(ceil(S*2*(cu-1)/t**2))
assert all(x>0 for x in dl)
# N(t) = (t^2 cosh(t/2)-2(cosh(t)-1))/t^4 has nonpositive
# coefficients beyond its constant 1/24, so is decreasing.
# D(t)=sinhc(t/2)^2 is increasing. S(t),A(t) likewise.
# Weight squared = 4 beta/(sinh(beta*(1-v^2/2))*sinhc(beta*v^2/2)).
# beta-grid width 1/40, v-grid width 1/256.
mins=[]
for i in range(70,480):
 lo=0
 for j in range(N):
  # t=beta*(1-v^2), interval rounded to scalar grid indices
  tl=G*i*(N*N-(j+1)**2)//(40*N*N)
  th=(G*(i+1)*(N*N-j*j)+40*N*N-1)//(40*N*N)
  num=nl[th];den=du[th] if num>=0 else dl[tl]
  c=floor(F(num*S,den))
  if c>=0:
   ai=(G*(i+1)*(2*N*N-j*j)+80*N*N-1)//(80*N*N)
   bi=(G*(i+1)*(j+1)**2+80*N*N-1)//(80*N*N)
   w=isqrt((i*S**4)//(10*su[ai]*au[bi]))
  else:
   ai=G*i*(2*N*N-(j+1)**2)//(80*N*N)
   bi=G*i*j*j//(80*N*N)
   w=isqrt(((i+1)*S**4)//(10*sl[ai]*al[bi]))+1
  lo+=c*w//S
 mins.append(F(lo,S*N)*F(2,5)) # sqrt(2)/pi>2/5; apply only to positive integrals
assert min(mins)>F(1,10**6), float(min(mins))
json.dump({'beta_range':['7/4','12'],'minimum_lower_exact':str(min(mins)),'minimum_approx':float(min(mins)),'cells':len(mins),'v_cells':N},open('M1_results.json','w'),indent=2)
print('PASS compact M1 lower bound',float(min(mins)))
