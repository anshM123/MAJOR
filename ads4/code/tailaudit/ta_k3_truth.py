"""
tailaudit/ta_k3_truth.py [DV]: how conservative is the binding k = 3 bound (TAIL_CERT §5: |H^(3)|/M <= 0.065 at L = 20)?
Same formula as tc_k3.py, |H^(3)|/M <= 0.2515 e^{-L} Lam^3/998^2 int_0^1 |psi3(p)||B(Lam p)| X3(p) dp / I_W(L), but with X3 replaced by
the actual fixed-twist |S^(7/2 + i Lam p)| at a = 998 (ta_power_check.Shat; regime-B path for Lam p > 519).  Not a bound: it shows
the size of the absolute-value k = 3 term with the true multiplier, i.e. how much of 0.065 is Theorem-P overestimate.
Output: cache/ta_k3_truth.txt
"""
import os, sys, json, math
from fractions import Fraction as Fr
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import ta_power_check as PC
from scipy.special import loggamma
ADS4 = PC.ADS4
nf = json.load(open(os.path.join(ADS4, 'prl', 'data', 'nf02_psi.json')))
K0 = Fr('27.48509817')
p3 = [float(Fr(c)/K0) for c in nf['psi3_exact']]
LAM = math.sqrt(998*999)
psi3 = lambda p: sum(c*p**i for i, c in enumerate(p3))
absB = lambda nu: math.sqrt(1/math.tanh(math.pi*nu)/(4*math.pi*nu))
NP = 400
ps = (np.arange(NP) + 0.5)/NP
I0 = 0.0; I1 = 0.0; Xmax = 0.0
for p in ps:
    nu = LAM*p
    v = abs(psi3(p))*absB(nu)/NP
    X = abs(PC.Shat(complex(3.5, nu), 998, 0, ft=True, Y1=(1.84*998 if nu > 0.52*998 else None)))
    Xmax = max(Xmax, X)
    I0 += v; I1 += v*X
W_I = json.load(open(os.path.join(ADS4, 'code', 'contour', 'cache', 'cb_w.json')))['cells']
L = 20.0
IW = sum(mn*(math.exp(-L*y0) - math.exp(-L*y1))/L for y0, y1, mn, mx_ in W_I)
pref = 0.2515*math.exp(-L)*LAM**3/998**2/IW
lines = [f'int_0^1 |psi3||B(Lam p)| dp ~ {I0:.4f} (doc: ~0.2);  with the actual fixed-twist |S^(7/2+i nu)| (max {Xmax:.3f}): {I1:.4f} '
         f'(certified X3-weighted integral: 329)',
         f'k = 3 ratio at L = 20 with the actual multiplier: {pref*I1:.3e}  (certified bound 0.065)']
for s in lines:
    print(s)
open(os.path.join(HERE, 'cache', 'ta_k3_truth.txt'), 'w').write('\n'.join(lines) + '\n')
