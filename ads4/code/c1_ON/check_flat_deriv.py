# checks: analytic u-derivatives of Pi0 vs numerical differentiation; large-p asymptotics vs OPE (T_mu nu) term
from flat import *
from mpmath import diff, nsum
mp.dps = 30
for n, u in [(0, mpf('0.7')), (0, mpf(5)), (1, mpf('2.3')), (3, mpf(40))]:
    f = lambda x: Pi0(n, x)
    print(n, u, 'd1', Pi0(n, u, 1), diff(f, u, 1), ' d2', Pi0(n, u, 2), diff(f, u, 2))
# OPE: eps = 8 p Pi - 1 -> -128 (h1-h2) (w^2 - p^2/3)/p^5
fR = lambda R: exp(-M0*R)/(4*pi*R)
f1 = lambda R: diff(fR, R, 1); f2 = lambda R: diff(fR, R, 2)
h1mh2 = nsum(lambda k: f2(k) - f1(k)/k, [1, inf])
print('h1-h2 =', h1mh2)
for n, u in [(0, mpf(400)), (0, mpf(1600)), (5, mpf(0)), (10, mpf(0)), (5, mpf(1000)), (10, mpf(4000))]:
    w = 2*pi*n; p = sqrt(u + w*w)
    eps = 8*p*Pi0(n, u) - 1
    pred = -128*h1mh2*(w*w - p*p/3)/p**5
    print(n, u, 'eps', eps, 'pred', pred, 'ratio', eps/pred)
