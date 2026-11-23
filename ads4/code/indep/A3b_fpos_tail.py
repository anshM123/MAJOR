"""F_{-1}(t) tail for the continuum C4 weights out to t = 22: |F| e^{3t/2} (-> const iff chi has a pole at -3i/2)."""
import os, sys, math
import mpmath as mp, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from A3_fpos import load_R, load_weight
import sympy as sp
mp.mp.dps = 100
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'A3b_fpos_tail.txt')
nu, R = load_R(); Rm1 = sp.lambdify(nu, R[-1], 'mpmath')
xs, _ = np.polynomial.legendre.leggauss(32)
xg, wg = [], []
for x0 in xs:
    x = mp.mpf(x0)
    for _ in range(6):
        x -= mp.legendre(32, x) / mp.diff(lambda u: mp.legendre(32, u), x)
    dp = mp.diff(lambda u: mp.legendre(32, u), x); xg.append(x); wg.append(2 / ((1 - x**2) * dp**2))
nodes, wts = [], []
for k in range(280):
    a, b = mp.mpf(k) / 4, mp.mpf(k + 1) / 4
    for x, w in zip(xg, wg):
        nodes.append((a + b) / 2 + (b - a) / 2 * x); wts.append((b - a) / 2 * w)
out = []
for piece in ('c4v', 'c4d'):
    chi, Q, B = load_weight(piece)
    base = [w * chi(v) for v, w in zip(nodes, wts)]
    Rp = [Rm1(v) for v in nodes]; Rn = [Rm1(-v) for v in nodes]
    for lab, fac in (('weight', lambda v: 1), ('control', lambda v: (v**2 + 4) / (v**2 + mp.mpf(9) / 4))):
        fs = [fac(v) for v in nodes]
        row = []
        for t in range(6, 23, 2):
            F = sum(bb * f * (a * mp.exp(-1j * v * t) + c * mp.exp(1j * v * t)) for v, bb, f, a, c in zip(nodes, base, fs, Rp, Rn))
            row.append((t, abs(F)))
        s = f'{piece} {lab:8s}: |F_-1(t)| e^(3t/2), t = 6..22 step 2: ' + ' '.join(mp.nstr(F * mp.exp(1.5 * t), 3) for t, F in row)
        rr = [(math.log(float(row[k][1])) - math.log(float(row[k + 2][1]))) / 4 for k in range(len(row) - 2)]
        s += '\n      mean decay rate over windows of 4: ' + ' '.join(f'{x:.3f}' for x in rr)
        print(s, flush=True); out.append(s)
open(OUT, 'w').write('\n'.join(out) + '\n')
