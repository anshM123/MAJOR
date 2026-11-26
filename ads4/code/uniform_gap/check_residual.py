#!/usr/bin/env python3
"""Exact arithmetic checks of the constants in the spectral-defect lemma.
The analytic bounds, including the singular-integral argument, are in RESEARCH_NOTE.md.
"""
from check_uniform import *
from math import factorial
p=c2[:]
for j in range(4):
 a=p[:]
 for _ in range(5-j):
  out=[];s=F(0)
  for x in a[:-1]:s+=x;out.append(s)
  assert a[-1]==-s
  a=out
 assert sum(abs(x) for x in a)<2**40
 p=[k*p[k] for k in range(1,len(p))]
# log F derivatives of P and cosh, and Gamma factors.
centers=[F(1,2),F(1,2),F(3,2),F(5,2),F(9,2)]
for j in (1,2,3):
 pol=sum(2*factorial(j-1)/c**j for c in centers)
 hyper=[0,F(11,7),F(11,7)**2,2*F(11,7)**3][j]
 gam=F(2) if j==1 else F(2*factorial(j-1),2**j)*(F(4,5)**j+F(4,5)**(j-1)/(j-1))
 assert pol+hyper+gam<100
assert 100**3+3*100**2+100<2**21
assert 26**13<2**62
# e^5>128 and e^6>400, exact series lower bounds.
assert sum(F(5)**k/factorial(k) for k in range(10))>128
assert sum(F(6)**k/factorial(k) for k in range(21))>400
assert F(12,399)+F(3,2)*(F(12,1600)+F(4,12))+F(36,121)<F(9,10)
# The combined defect norm estimate is bounded by 2^212(e^-5000+1002e^-1000).
assert 212+10+1-7*200<-1000
# Part L, M1 compact floor from independent monotone-table check.
M=json.load(open(Path(__file__).with_name('M1_results.json')))
assert F(M['minimum_lower_exact'])>F(1,10**6)
assert F(3,2**1000)<F(1,10**6)
# Light bound constants.
L=exact_light()
assert 16*F(L['slope'])<3023
assert F(L['E_constant_upper'])+F(1,2**500)<-5307
assert F(L['K_min'])-F(1,2**500)>F('0.99995')
print('PASS: residual derivative and endpoint constants, compact M1 check, analytic tail M1 lower bound, logarithmic light-bound constants')
