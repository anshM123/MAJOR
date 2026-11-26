#!/usr/bin/env python3
"""Exact contribution of a rational-moment seed error to the direct prefix.
Raw moment pole at nu^2=-a/4 is an imaginary spectral atom. The whole
coefficient sequence is a polynomial in k (apart from a trivial k=0 check).
No high-order moment error amplification is allowed for these correlated errors.
"""
from fractions import Fraction as F
import sympy as s
import json
k,z=s.symbols('k z');R=F(99,100);C=F(9,10);hmax=F(1,998**2)
polys=[]
for a,l in [(1,0),(9,1),(49,3)]:
  P=s.cancel(s.legendre(l,(1+z)/(1-z))*(1-z)**l).expand();co=s.Poly(P,z)
  # [z^(m-1)] F_2: expression for m>=1, binomial polynomial vanishes
  # in each negative-index case (degree l+1>=needed shifts).
  m=s.symbols('m');g=sum(co.nth(j)*s.prod(m-1-j+i for i in range(1,l+2))/s.factorial(l+1) for j in range(l+1));f=s.cancel(g/m**2)
  def Z(j):return s.cancel(12*f.subs(m,j+2)-(j+1)*f.subs(m,j+1)+(j+3)*f.subs(m,j+3))
  D=(k+1)*(k+2)*(k+3)/6
  A=s.cancel(D**2*(Z(k)-k*Z(k-1)/(k+3))/(k+3));A=s.Poly(A,k)
  assert A.eval(0)==Z(0)/3
  # independent exact recurrence at low indices (entire atom action)
  fs=[F(0),F(1)]
  for n in range(1,30):fs.append(((2*n*n+F(a-1,4))*fs[-1]-(n-1)**2*fs[-2])/(n+1)**2)
  zz=[12*fs[j+2]-(j+1)*fs[j+1]+(j+3)*fs[j+3] for j in range(25)]
  from math import comb
  for j in range(25):assert F(A.eval(j))==F(comb(j+3,3)**2,j+3)*(zz[j]-F(j,j+3)*(zz[j-1] if j else 0))
  # A and dressed coefficients polynomial. Sum absolute coefficient l1
  # against c^k using exact power-geometric moments; physical U_k <= c^k.
  for b in (0,4):
    q=b*(b+1);Dpol=s.Poly(((k+1)**2*A.as_expr().subs(k,k+1)+(-k*k+(8-4*b)*k+9-2*b*(b-1))*A.as_expr())/2,k)
    def bound(P):
      c=s.Rational(C.numerator,C.denominator);gs=1/(1-z);tot=0
      for j in range(P.degree()+1):
        if j:gs=s.cancel(z*s.diff(gs,z))
        tot+=abs(P.nth(j))*gs.subs(z,c)
      return F(tot)
    # full physical normalization w^3<=1 and tau<=h<=1/998^2
    norm=bound(A)+hmax*bound(Dpol)
    polys.append(dict(a=a,b=b,A=str(A.as_expr()),dress=str(Dpol.as_expr()),action_bound=str(norm),approx=float(norm)))
    print(a,b,'A',A.as_expr(),'bound',float(norm))
# R partial fractions. g = -(1/Lambda) R B theta, F_a is half-line base.
coeff={1:F(-96),9:F(5184,5),49:F(-4704,5)}
# Actual inventory bound S_cut < 7.1e-24, proved separately by exact integer sums.
# |Delta F_a| <= S_cut/a. (The log term is not covered by this bound.)
for b in (0,4):
  err=F('7.1e-24')/998*sum(abs(coeff[r['a']])*F(r['action_bound'])/r['a'] for r in polys if r['b']==b)
  assert err<F('1e-10');print('cutoff pole contribution',b,float(err))
json.dump(polys,open('pole_bound_full_results.json','w'),indent=2)
