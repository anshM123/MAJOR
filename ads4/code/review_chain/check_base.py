#!/usr/bin/env python3
"""Referee check: base_bound (M0,M1 envelopes for |X(alpha)|,|X(alpha+1)|) against exact Arb ladder sums
on physical states at the edges of the L-cells, at real and complex (top-line, negative Re alpha) frequencies.
usage: check_base.py CELL   (CELL in 26, 29)"""
import sys
from fractions import Fraction as F
from flint import arb, acb, ctx
sys.path.insert(0, '/home/user/MAJOR/ads4/code/sliver_ext')
from base_bounds import base_bound
ctx.prec = 9000

def solve_n(a, b, d_target):
    s = 2*b+1; lo, hi = 1, 10**7
    while hi-lo > 1:
        m = (lo+hi)//2; N = m+a
        d = F(a*(N+m+s), N*(N+s))
        if d <= d_target: hi = m
        else: lo = m
    return hi

def Xexact(n, a, b, alpha):
    s = 2*b+1; N = n+a; D = N*(N+s)
    c = arb(n*(n+s))/D; h = arb(1)/D; d = 1-c
    tot = acb(0); term = acb(1); U = arb(1)
    for k in range(n+1):
        tot += term*U
        term = term*(alpha+k)/(k+1)
        U = U*(c-h*k*(k+s))/(1-h*k*(k+s))
    return abs(acb(d)**alpha*tot), float(F(a*(N+n+s), D))

if __name__ != "__main__": raise SystemExit
cells = {'26': (F(68779, 1250000), F(46421, 625000)), '29': (F('0.0407622')-F('0.0015'), F('0.0407622')+F('0.0015'))}
dl, du = cells[sys.argv[1] if len(sys.argv) > 1 else '26']
worst = 0
PS = [(3, ['0', '10.25']), (4, ['0', '20.25'])]
if len(sys.argv) > 2: PS = [x for x in PS if x[0] == int(sys.argv[2])]
for p, Ys in PS:
    for b in (0, 4):
        for Y in Ys:
            for t in ['0.25', '1', '160', '600', '950', '998.5']:
                t_ = F(t); y_ = F(Y)
                M0, M1 = base_bound(p, b, dl, du, t_, t_, y_, y_)
                for a in (998, 2500):
                    for dt in (dl, du):
                        n = solve_n(a, b, dt)
                        if not (dl <= F(a*(2*n+a+2*b+1), (n+a)*(n+a+2*b+1)) <= du):
                            n -= 1
                        q=F(2*p-1, 2)-y_; al = acb(arb(q.numerator)/q.denominator, arb(t_.numerator)/t_.denominator)
                        x0, d = Xexact(n, a, b, al)
                        x1, _ = Xexact(n, a, b, al+1)
                        assert float(x0.rad()) < 1e-6*max(1,float(x0.mid())) and float(x1.rad()) < 1e-6*max(1,float(x1.mid())), 'precision'
                        r0 = float(x0.upper())/float(M0); r1 = float(x1.upper())/float(M1)
                        worst = max(worst, r0, r1)
                        print(p, b, Y, t, a, n, round(d, 6), 'X0', float(x0.mid()), 'M0', float(M0), 'r0 %.3g' % r0, 'X1', float(x1.mid()), 'M1', float(M1), 'r1 %.3g' % r1, flush=True)
print('worst exact/bound', worst, 'PASS' if worst < 1 else 'FAIL')
