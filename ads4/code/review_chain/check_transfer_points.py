#!/usr/bin/env python3
"""Referee check: exact signed transfer coefficients P,Q (J0=1500, dressed for p=3) at physical states and
complex frequencies on the [2.6,2.9] contour/residual, compared with the certified bounds the assembly uses
(Extend.get = real norming + Hadamard ellipse extension, or the absolute majorant)."""
import sys, json, random
import mpmath as mp
from fractions import Fraction as F
sys.path.insert(0, '/home/user/MAJOR/ads4/code/sliver_next')
import os; os.chdir('/home/user/MAJOR/ads4/code/sliver_next')
from ellipse_extend import Extend
from check_norm import parameters
mp.mp.dps = 60
J = 1500

def PQ(d, h, nu, p, b):
    s = 2*b+1; q = b*(b+1)
    alpha = mp.mpf(2*p-1)/2+1j*nu
    def r(j):  # c_{j+1}/c_j
        a0 = mp.mpf(1)/2-1j*nu+j
        return a0*a0/((1-2j*nu+j)*(j+1))
    out = []
    for Z0, Zm in [(1, 0), (0, 1)]:
        Zp, Zc = mp.mpc(Zm), mp.mpc(Z0); rprev = mp.mpf(1); tot = 0
        for j in range(J):
            beta = alpha-j
            if p == 3:
                A0 = mp.mpf(9)/2-b*(b-1)-beta*(beta+8-4*b)/2; A1 = beta*(beta+9-4*b)/2
                tot += ((1+q*h+h*A0)*Zc+h*A1*rprev*Zp)/(1+q*h)
            else:
                tot += Zc
            Dj = 1-h*beta*(beta-s); rj = r(j)
            Zn = (rj*(d+h*beta*(s-2*beta-1))*Zc+rj*rprev*h*beta*(beta+1)*Zp)/Dj
            Zp, Zc, rprev = Zc, Zn, rj
        out.append(tot)
    return abs(out[0]), abs(out[1])/d   # Q = Q_unscaled/d

def phys(a, b, dt):
    s = 2*b+1; lo, hi = 1, 10**8
    while hi-lo > 1:
        m = (lo+hi)//2; N = m+a
        if F(a*(N+m+s), N*(N+s)) <= dt: hi = m
        else: lo = m
    n = hi; N = n+a; D = N*(N+s)
    return mp.mpf(a*(N+n+s))/D, mp.mpf(1)/D, n

rows = json.load(open('results_2.6_2.9.json'))
d0, dd, *_ = parameters(rows[0]['config']); dl, du = d0-dd, d0+dd
E = Extend(rows, dl, du)
random.seed(11)
worst = 0
pts = [(3, F('500'), F('10.25')), (3, F('900'), F('10.25')), (3, F('998.5'), F('5')), (3, F('1000'), F(0)), (3, F('300'), F(0)),
       (4, F('998.5'), F('12')), (4, F('600'), F('20.25')), (4, F('50'), F('20.25')), (3, F('5'), F('10.25'))]
for p, t, y in pts:
    for b in (0, 4):
        Pb, Qb = E.get(p, b, t, t, y, y)
        for a, dt in [(998, dl), (998, du), (1500, (dl+du)/2), (40000, du)]:
            d, h, n = phys(a, b, dt)
            if not (dl <= F(mp.nstr(d, 30)) <= du): continue
            nu = mp.mpc(float(t), float(y)) if False else mp.mpf(t.numerator)/t.denominator+1j*mp.mpf(y.numerator)/y.denominator
            P, Q = PQ(d, h, nu, p, b)
            rp, rq = P/(mp.mpf(Pb.numerator)/Pb.denominator), Q/(mp.mpf(Qb.numerator)/Qb.denominator)
            worst = max(worst, rp, rq)
            print(p, b, float(t), float(y), a, n, 'd=%.6f' % float(d), '|P|=%.4g bound %.4g' % (P, float(Pb)), '|Q|=%.4g bound %.4g' % (Q, float(Qb)), 'ratios %.3g %.3g' % (rp, rq), flush=True)
print('worst exact/bound', mp.nstr(worst, 5), 'PASS' if worst < 1 else 'FAIL')
