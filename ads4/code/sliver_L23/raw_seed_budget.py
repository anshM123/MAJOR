#!/usr/bin/env python3
"""Conservative raw real-frequency kernel bound for 0.05<=c<=0.56.
A seed error in K0 gives an actual Gaussian moment sequence and three
low imaginary-frequency pole atoms, not independent moment balls.
"""
from fractions import Fraction as F
import sympy as s
z,k=s.symbols('z k');R=F(99,100);r=F(14,25)/R
# f_m = [z^(m-1)]F_2/m^2, so |f_m| <= M R^-(m-1)/m^2,
# M=20000 exp(11 nu/7). Then |Z_k|<=5 M R^-(k+2),
# |A_k|<=10 M R^-(k+2)(k+3)^5, and
# |B_k|<=20 M R^-(k+3)(k+4)^7 (both channels).
def geomsum(P):
    c=s.Rational(r.numerator,r.denominator);p=s.Poly(P,k);g=1/(1-z);v=0
    for j in range(p.degree()+1):
        if j:g=s.cancel(z*s.diff(g,z))
        v+=abs(p.nth(j))*g.subs(z,c)
    return F(v)
C=20000*(10/R**2*geomsum((k+3)**5)+20/R**3/998**2*geomsum((k+4)**7))
assert C<10**14
# P(nu)=(1+nu^2)sum y_i nu^(2i), |y_i|<1, so |P|<=3(1+nu)^6.
# |R_t| <= 9216/45+55296/245 < 500.
assert F(9216,45)+F(55296,245)<500
# e^(-nu^2/2+11nu/7)<= e^(121/49)e^-nu^2/4 < 16 e^-nu^2/4.
from check_lowfreq import exp_u
assert exp_u(F(121,49))<16
# (1+nu)^6<=32(1+nu^6), integral_0^inf e^-nu^2/4 = sqrt(pi)<2,
# integral_0^inf nu^6 e^-nu^2/4 =120 sqrt(pi)<240.
Q=3*32*16*(2+240);assert Q<10**6
# If K0 error <=1e-60, G0>1, the induced Gaussian polynomial rational-weight error is below:
E=F('1e-60')*C*500*Q/998
# Pole seed correction of K0-induced particular solution <=3*32*(2+240)*1e-60/a;
# actual LF0 seeds accurate to 1e-60. Use exact pole bounds from pole_bound.py.
import json
rr=json.load(open('pole_bound_results.json'));coef={9:F(-9216,5),49:F(55296,5)}
for b in (0,4):
    e=E+sum(abs(coef[x['a']])*F(x['action_bound'])*F('1e-60')*(1+F(3*32*(2+240),x['a']))/998 for x in rr if x['b']==b and x['a'] in coef)
    assert e<F('1e-38');print('b',b,'raw kernel constant',float(C),'correlated log seed action error',float(e))
# Dropping theta on the log component only: on [0,64], 1-theta < exp(-10000).
# On nu>=64, log(1+nu^2)<=nu^2 and
# 3*(1+nu)^6*nu^2 <=192*nu^8 <= exp(nu^2/12), by logarithmic derivative
# (8/nu <= nu/6) and the inequality at nu=64.
from base_bounds import logu,logl
assert logu(192)+8*logu(64)<F(64**2,12)
assert F(11,7)<F(64,6)
# Weighted Gaussian e^(-nu^2/2+11nu/7) <= e^(-nu^2/3) there.
# Thus tail integral <= int_64^inf e^(-nu^2/4)<=e^-1024/32.
# On [0,64], additional nu^2<=4096, and use Q above.
from resid_budget import pow2exp
cut=C*500/998*(4096*Q*pow2exp(-10000)+pow2exp(-1024)/32)
assert cut<F('1e-400');print('Log cutoff total action error <1e-400')
