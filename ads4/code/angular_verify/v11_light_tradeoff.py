"""v11: how small must the light remainder be?  Add ALL light states (twist 2..g, both ladders, J <= 4N+8; the light
[TT]-type double twists) to the LP as columns whose endpoint cost is kappa*W instead of W, i.e. the light share of
the endpoint budget is allowed to be E_light <= 3/kappa (A2' says E_light -> 0 as c_T -> infinity).
Prints |c|^2_max (and tau g^3) versus kappa.   usage: python3 v11_light_tradeoff.py g N"""
import sys, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
import v3_lp_hp as V
from common import allowed
g = int(sys.argv[1]); N = int(sys.argv[2])
M = V.Mmat(N)
MT = V.arb_mat([[V.arb(M[n][k].numerator) / M[n][k].denominator for n in range(N + 1)] for k in range(N + 1)])
S = V.candidates(g)
S = [s for s in S if s[1] <= 3 * g + 10 and float(s[2]) - s[1] <= 30 * g]
Lt = [Fr(2), Fr(5, 2), Fr(3), Fr(4), Fr(6)] + [Fr(x).limit_denominator(16) for x in np.linspace(7, g - 0.5, 12)]
SL = [(b, J, J + t) for t in Lt for J in range(0, 4 * N + 9) for b in (0, 4) if allowed(b, J)]
B = V.basis_rows(S, N, MT); BL = V.basis_rows(SL, N, MT)
a_eq = np.array([float(M[n][2]) for n in range(N + 1)])
for kappa in [1e30, 1e20, 1e15, 1e12, 1e9, 1e6, 1e3, 1.0]:
    R = np.vstack([B, BL / kappa])
    A = np.hstack([R, -np.ones((R.shape[0], 1))]); c = np.zeros(N + 2); c[-1] = 1
    res = linprog(c, A_ub=A, b_ub=np.zeros(R.shape[0]), A_eq=np.append(a_eq, 0)[None, :], b_eq=[100.0],
                  bounds=[(None, None)] * (N + 1) + [(0, None)], method='highs')
    if res.status != 0:
        print(f"  kappa={kappa:.0e}: LP status {res.status} ({res.message[:40]})"); continue
    c2 = 3 * res.x[-1]
    print(f"  g={g} N={N} kappa={kappa:.0e} (E_light <= {3/kappa:.0e}): |c|^2 <= {c2:.4e}, tau <= {4*np.sqrt(c2):.4e}, tau g^3 = {4*np.sqrt(c2)*g**3:.1f}")
