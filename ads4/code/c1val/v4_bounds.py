#!/usr/bin/env python3
"""
c1val / v4: Result D and Result Q as functions of c1, from the exact closed forms of C1_COMPUTE.md sec. 4
(constants copied from code/c1/cache/c1_quadratic.txt; nothing re-derived here).
  D:  sqrt(t4^2 + aT^2) <= R0 (x + sqrt(1 + x^2)),   x = xi |c1|
  Q:  g4~ Lam^6 / 8 pi G <= 2 L_AdS + q c1^2
Also: the |c1| at which each bound moves by 1 % / 10 %, and the values at the natural size and at the C1_VALUE premises.
Writes cache/v4_bounds.txt.
"""
import os, math

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'v4_bounds.txt')
DATA = {300: dict(R0=3.13206e-6, xi=2.16006e-3, twoL=72.04059, q=6.81453e-5),
        1000: dict(R0=2.89126e-8, xi=1.67026e-4, twoL=92.81351, q=5.44276e-7)}
lines = []
def log(s=''):
    print(s); lines.append(s)

def Dratio(x):
    return x + math.sqrt(1 + x*x)

log('v4: exact c1-dependence of D and Q (constants from code/c1/cache/c1_quadratic.txt)')
log('| Delta_g | |c1| | x | D bound sqrt(t4^2+aT^2) <= | D ratio | Q upper end | Q shift (rel.) |')
log('|---|---|---|---|---|---|---|')
for dg, p in DATA.items():
    for c1 in [0, 1e-3, 1e-2, 0.1, 1, 10, 100]:
        x = p['xi']*c1
        Rb = p['R0']*Dratio(x)
        Qu = p['twoL'] + p['q']*c1*c1
        log(f'| {dg} | {c1:g} | {x:.3e} | {Rb:.5e} | {Dratio(x):.6f} | {Qu:.5f} | {p["q"]*c1*c1/p["twoL"]:.2e} |')
log('')
log('|c1| at which a bound moves by a given relative amount:')
for dg, p in DATA.items():
    for rel in [1e-3, 1e-2, 1e-1]:
        # D: x + sqrt(1+x^2) = 1 + rel  ->  x = ((1+rel)^2 - 1) / (2 (1+rel))
        r = 1 + rel
        xD = (r*r - 1)/(2*r)
        cD = xD/p['xi']
        cQ = math.sqrt(rel*p['twoL']/p['q'])
        log(f'  Delta_g = {dg}: {rel*100:g} %  ->  D: |c1| = {cD:.4g}    Q: |c1| = {cQ:.4g}')
with open(OUT, 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
