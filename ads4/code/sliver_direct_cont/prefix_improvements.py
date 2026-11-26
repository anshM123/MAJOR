#!/usr/bin/env python3
"""Exact bounds for the original and conformally transformed physical prefix tails.
Uses the spectral mass and Gaussian envelopes retained in README.md.
"""
from fractions import Fraction as F
from base_bounds import logl
from resid_budget import gaussian_mass
import math
H=F(1,998**2)
V=1040
B=F(11,7)
R=F(199,200)
assert V*B<2358*logl(F(2))
M=2*(1-R)**-4*(120000001*2**2358+F(2,998)*gaussian_mass(V,B))
assert M<2**2417

def tail(c,K):
 t=F(c)/R;assert 0<t<1
 return 2**2417*t**(K+1)*((1+20*H)/(1-t)+2*H*((K+5)**2/(1-t)+2*(K+5)*t/(1-t)**2+t*(1+t)/(1-t)**3))
TIGHT_ALLOW=F(3,10**6)
assert tail(F('.258'),1251)<TIGHT_ALLOW
# Conformal map w=(1-sqrt(1-z))/(1+sqrt(1-z)). On |w|<=.99,
# (1-z)^-p=((1+w)/(1-w))^(2p) and beta=4 atanh(sqrt(w)).
Rc=F(99,100);wc=F(13,25);Bc=F(22,7)
assert (F(6,19))**2<F(1,10) # w(.9)<.52
assert 4*F(99,100)*(1-Rc)/(1+Rc)**2>F(99,10000)
assert Bc/F(99,10000)<320
assert V*Bc<4716*logl(F(2))
Mc=320*199**8*(120000001*2**4716+F(2,998)*gaussian_mass(V,Bc))
assert Mc<2**4813
# w' = 1/[sqrt(1-c)(1+sqrt(1-c))^2]; w''=1/[2 t^3(1+t)^2]+1/[t^2(1+t)^3].
assert 1/(F(6,19)*(1+F(6,19))**2)<4
assert 1/(2*F(6,19)**3)+1/F(6,19)**2<26
assert 16+26+20+F(15,2)<70 # used only as a loose auxiliary check
# More precisely f''+5f'+7.5f <= (16j(j-1)+46j*wc+7.5*wc^2)*wc^(j-2)
# <= 50*(j+1)^2*wc^(j-2), termwise.
assert 46*wc<100 and F(15,2)*wc**2<50
assert (1+20*H)*wc**2+50*H<3

def conformal_tail(J):
 t=wc/Rc
 return 3*2**4813*wc**-2*t**(J+1)*((J+2)**2/(1-t)+2*(J+2)*t/(1-t)**2+t*(1+t)/(1-t)**3)
assert conformal_tail(5333)<F(1,2**110)
if __name__=='__main__':
 for name,x in [('direct_tail_at_.258',tail(F('.258'),1251)),('conformal_tail_at_.9',conformal_tail(5333))]:
  print(name, float(x),'log2',math.log2(x.numerator)-math.log2(x.denominator))
 print('PASS_EXACT_TAIL_BOUNDS')
