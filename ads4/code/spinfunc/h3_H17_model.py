"""
h3 (M5'b): the |ghat3|^2 coefficient k_AdS = k(0)(1 + dk) for route A' in the SAME heuristic model as omegau/nf5_light.py (b) [M]:
each light-side |ghat3|^2 polynomial of degree D in q = p^2 receives AdS corrections with coefficients |a| = D x 36 per power lowered
(CHMRS-type subleading terms, cf. -C_{2,nu}|_{g3} = nu^2 - (2Df - (d-4)/2)(2Df + (d-4)/2)); dk ~ sum_k |lev_k| D_k 36/L^2,
lev_k = (coefficient of q^D)/pi x int_0^1 psi_k p^{2D-2} dp.
Route A': |ghat3|^2 only in the stripped B2 (-B2 ∋ 2 pi G |g3|^2 q^3) and B3 (-B3 ∋ -2 pi G |g3|^2 q^2), weights psi_2, psi_3.
Old Omega^{NF,SC} (B' uplift): B'_2 (q^3), B'_3 (q^2), B'_4 (2 q) with weights phi_k (reproduces nf5's 0.042 at Dg = 100, n = 1).
Effect: A, B_AdS -> /(1 + dk); bounds -> (1 + dk)^{-1/2}; signs of the true corrections unknown -> quoted as +-.
Usage: python3 h3_H17_model.py -> cache/h3_H17_model.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
import numpy as np, sympy as sp
from nf_common import NF
F = NF(); p = sp.Symbol('p')
PS = {2: sum(sp.Rational(c.numerator, c.denominator)*p**i for i, c in enumerate(F.s2)),
      3: sum(sp.Rational(c.numerator, c.denominator)*p**i for i, c in enumerate(F.s3))}
PH = {k: sum(sp.Rational(c.numerator, c.denominator)*p**i for i, c in enumerate(F.PHx[k])) for k in range(2, 8)}
mu = lambda f, j: float(sp.integrate(sp.expand(f*p**j), (p, 0, 1)))
# (coefficient of q^D in -B|_{|g3|^2}/(pi G |g3|^2), degree D)
routeA = {2: (2, 3), 3: (-2, 2)}                         # stripped B2, B3
oldBp = {2: (2, 3), 3: (2, 2), 4: (4, 1)}                # unstripped B'_2, B'_3, B'_4 (FLAT_UNSTRIPPED table)
def dk_of(table, W, L):
    return sum(abs(c*mu(W[k], 2*D - 2))*D for k, (c, D) in table.items())*36/L**2
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)
pr("# h3: H17 model (nf5 heuristic, |a| = D x 36): relative k_AdS - k(0) at n = 1, and the induced bound factor (1 -+ dk)^(-1/2)")
for Dg in (100, 300, 1000):
    L = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    dA = dk_of(routeA, PS, L); dB = dk_of(oldBp, PH, L)
    pr(f"  Delta_g = {Dg}: route A' dk ~ {dA:.3g}  (old B'-uplift, same model: {dB:.3g});  bounds x [{(1 + dA)**-0.5:.4f}, {(1 - dA)**-0.5:.4f}];"
       f"  A, B_AdS x [{1/(1 + dA):.4f}, {1/(1 - dA):.4f}]")
open(os.path.join(HERE, 'cache', 'h3_H17_model.txt'), 'w').write("\n".join(out) + "\n")
