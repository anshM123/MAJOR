"""
Referee check (d): the bound formula  sqrt(t4^2+aT^2) <= 240 sqrt(A log(Lam/n) + B_AdS) (n/Lam)^4  and the exact L_AdS form used by e6/e9,
plus the dictionary arithmetic (240 = sqrt(57600); ECG cross-check t4 = -1260 mu).
Independent pieces: L_AdS integral done here by mpmath on the exact NF02/M02 psi2 polynomials; A, B_AdS re-derived from psi2 analytically:
  16 int_0^1 psi2(p)/(p^2+e^2) dp  =  A log(1/e) + B0 + O(e),  A = 16 psi2'(0),  with e = 3n/(2 Lam)
  => A log(Lam/n) + [B0 + A log(2/3)] = A log(Lam/n) + B_AdS.
Output: cache/rf2_bound.txt
"""
import os, sys, json
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
mp.mp.dps = 30
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

from nf_common import NF
F = NF()
# recover NF02 psi2 as a polynomial by exact fit on many points (it is a polynomial of moderate degree)
ps = np.linspace(0, 1, 401)
vals = np.array([float(F.psi2(p)) for p in ps])
deg = None
for d in range(3, 30):
    c = np.polynomial.polynomial.polyfit(ps, vals, d)
    if np.max(np.abs(np.polynomial.polynomial.polyval(ps, c) - vals)) < 1e-9*np.max(np.abs(vals)):
        deg = d; break
pr(f"# rf2: NF02 psi2 is a polynomial of degree {deg} (fit residual < 1e-9 rel); psi2(0) = {vals[0]:.3e}")
psi2_nf = lambda p: mp.mpf(float(F.psi2(float(p))))
d0 = float(mp.diff(lambda p: float(F.psi2(float(p))), mp.mpf('1e-6')))
Jd = json.load(open(os.path.join(HERE, '..', 'robust', 'functionals_certified.json')))
cM = [Fr(x) for x in Jd['M02']['psi2_exact']]; k0M = Fr(10884217, 396000)
psi2_m02 = lambda p: sum(mp.mpf(ci.numerator)/ci.denominator*p**i for i, ci in enumerate(cM))

def LAdS(psi2, L):
    e2 = mp.mpf(9)/(4*L*L)
    return mp.quad(lambda t: psi2(t)/(t*t + e2), [0, mp.sqrt(e2), 0.1, 1])

def AB(psi2, k0):
    # A from psi2'(0); B0 from 16 int (psi2 - psi2'(0) p)/p^2 + 16 psi2'(0) int_0^1 p/(p^2) regularised -> use large-L fit instead
    Ls = [mp.mpf(10)**j for j in (5, 6, 7)]
    ys = [16*LAdS(psi2, L)/k0 for L in Ls]
    A = (ys[2] - ys[1])/mp.log(10)
    B = ys[2] - A*mp.log(Ls[2])
    return A, B

for name, psi2, k0, Aclaim, Bclaim in (('NF02 (D-safe)', psi2_nf, 1, 3022.3377, -6552.6620),
                                       ('M02 (D-strong)', psi2_m02, mp.mpf(k0M.numerator)/k0M.denominator, 51.7366, -86.2981)):
    A, B = AB(psi2, k0)
    pr(f"## {name}: A = {float(A):.4f} (claim {Aclaim}), B_AdS = {float(B):.4f} (claim {Bclaim}); B_AdS - A log(2/3) = flat B = {float(B - A*mp.log(mp.mpf(2)/3)):.4f}")
    for Dg, n in ((100, 1.0), (300, 1.0), (1000, 1.0), (100, 1.08), (300, 1.05), (1000, 1.01)):
        Lam = mp.sqrt((Dg - 2)*(Dg - 1))
        X = n**8*16*LAdS(psi2, Lam/n)/k0
        exact = 240*mp.sqrt(X)/Lam**4
        lead = 240*mp.sqrt(A*mp.log(Lam/n) + B)*(n/Lam)**4
        pr(f"   Dg={Dg} n={n}: |g3|^2 Lam^8 <= {float(X):.5g};  sqrt(t4^2+aT^2) <= {float(exact):.4g} (exact L_AdS)  vs  {float(lead):.4g} (A log + B_AdS form)")
    # sensitivity: d log(bound)/d log n
    Dg = 100; Lam = mp.sqrt(98*99)
    f = lambda n: 240*mp.sqrt(n**8*16*LAdS(psi2, Lam/n)/k0)/Lam**4
    pr(f"   sensitivity at Dg=100: bound(n=1.2)/bound(1) = {float(f(1.2)/f(1.0)):.3f}, bound(1.5)/bound(1) = {float(f(1.5)/f(1.0)):.3f}, bound(2.0)/bound(1) = {float(f(2.0)/f(1.0)):.3f}")

pr("## dictionary arithmetic")
pr(f"   sqrt(57600) = {mp.sqrt(57600)};  240 = 2(d-1)d(d+1)(d+2) at d=3: {2*2*3*4*5}")
# ECG: L = -(mu L^4/8) P, P ~ 12 I2 + I1 = 7 I1 on the 3-pt vertex (I2 = I1/2); CLPS: -(alpha3/6) I1  => alpha3 = 6*7/8 mu = 21/4 mu
pr(f"   ECG: alpha3 = 6*7/8 mu = {Fr(42,8)} mu;  t4 = -240 alpha3 = {-240*Fr(42,8)} mu  (BCR 1802.00018: t4 = -1260 mu at leading order)")
# QTG d=4: (7/8) mu L^4 Z5', Z5' = I1 + ...  => -(alpha3/6) = 7/8 mu => alpha3 = -21/4 mu, t4 = -720 alpha3 = 3780 mu
pr(f"   QTG d=4: alpha3 = -6*7/8 mu = {-Fr(42,8)} mu;  t4 = -720 alpha3 = {720*Fr(42,8)} mu  (MPS 1004.2055: 3780 mu)")
open(os.path.join(HERE, 'cache', 'rf2_bound.txt'), 'w').write("\n".join(out) + "\n")
