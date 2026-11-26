#!/usr/bin/env python3
"""Referee check: B_box / B_bound / inf_BH / B_inf envelopes against mpmath on complex nu=t+iy,
plus the identity |Gamma(2i nu)Gamma(b0)/Gamma(a0)^2/Gamma(1/2+i nu)^2| = |coth(pi nu)|/(2 pi)."""
import sys, random
import mpmath as mp
from fractions import Fraction as F
sys.path.insert(0, '/home/user/MAJOR/ads4/code/sliver_ext')
from contour_direct import B_box, inf_BH
from assemble_top import B_bound, B_inf
mp.mp.dps = 40
random.seed(7)

def B(nu): return mp.gamma(2j*nu)/mp.gamma(0.5+1j*nu)**2
def H(nu, d):
    a0 = 0.5-1j*nu; b0 = 1-2j*nu
    return mp.hyp2f1(a0, a0, b0, d)

worst_B = worst_BH = 0; cnt = 0
pts = []
for t in ['0.25', '0.3', '0.5', '1', '3', '10', '50', '159', '400', '998.5']:
    for y in ['0', '0.25', '0.49', '0.5', '0.51', '0.99', '1', '1.01', '2.5', '5', '6.25', '10.25', '11.25', '19.25', '20.25']:
        pts.append((F(t), F(y)))
for _ in range(150):
    pts.append((F(random.randint(25, 99850), 100), F(random.randint(0, 2100), 100)))
for t, y in pts:
    nu = mp.mpc(mp.mpf(t.numerator)/t.denominator, mp.mpf(y.numerator)/y.denominator)
    bb = abs(B(nu)); env = F(B_box(t, t, y, y)); envf = mp.mpf(env.numerator)/env.denominator
    worst_B = max(worst_B, bb/envf)
    if y >= F(1, 2) or True:
        if y > 0:
            envt = F(B_bound(t, t, y)); worst_B = max(worst_B, bb/(mp.mpf(envt.numerator)/envt.denominator))
    for d in ['0.026', '0.0550232', '0.0742736']:
        dv = mp.mpf(d); du = F(d)
        bh = abs(B(nu)*H(nu, dv))
        e2 = F(inf_BH(du, y, y)); r = bh/(mp.mpf(e2.numerator)/e2.denominator)
        worst_BH = max(worst_BH, r)
        if y >= F(1, 2):
            e3 = F(B_inf(du, y)); worst_BH = max(worst_BH, bh/(mp.mpf(e3.numerator)/e3.denominator))
    cnt += 1
print('points', cnt, 'worst |B|/env', mp.nstr(worst_B, 6), 'worst |BH|/env', mp.nstr(worst_BH, 6))
# interval version: B_box over [t0,t1]x[y0,y1] vs sampled max
worst_i = 0
for (t0, t1, y0, y1) in [(F(1,4), F(1,4), F(0), F(1,4)), (F(1,4), F(1,4), F(3,4), F(1)), (F(1,4), F(1,4), F(19,2), F(39,4)),
                         (F(1,4), F(1,2), F(41,4), F(41,4)), (F(1), F(10), F(81,4), F(81,4)), (F(1997,2), F(1997,2), F(10), F(41,4))]:
    env = F(B_box(t0, t1, y0, y1)); m = 0
    for i in range(9):
        for k in range(9):
            t = t0+(t1-t0)*i/8; y = y0+(y1-y0)*k/8
            nu = mp.mpc(mp.mpf(t.numerator)/t.denominator, mp.mpf(y.numerator)/y.denominator)
            m = max(m, abs(B(nu)))
    worst_i = max(worst_i, m/(mp.mpf(env.numerator)/env.denominator))
    print('box', [float(x) for x in (t0, t1, y0, y1)], 'max|B|', mp.nstr(m, 5), 'env', float(env))
print('worst interval ratio', mp.nstr(worst_i, 5))
# identity coth/(2pi)
for nu in [mp.mpc(0.3, 2.2), mp.mpc(5, 7.5), mp.mpc(900, 10.25)]:
    a0 = 0.5-1j*nu; b0 = 1-2j*nu
    lhs = abs(mp.gamma(2j*nu)*mp.gamma(b0)/mp.gamma(a0)**2/mp.gamma(0.5+1j*nu)**2)
    print('euler prefactor', mp.nstr(lhs, 12), mp.nstr(abs(mp.coth(mp.pi*nu))/(2*mp.pi), 12))
print('PASS' if worst_B < 1 and worst_BH < 1 and worst_i < 1 else 'FAIL')
