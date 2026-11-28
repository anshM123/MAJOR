"""v6: size of the optimal all-k dual functional on LIGHT states (twist < g), which sets the coefficient of the
A2' light remainder in the all-k bound:  |c|^2 <= 3 lam + sum_light rho (sigma f - lam W)_+.
For the paper's k=2 rule the analogous coefficient is |a_2|^2/|a_0|^2 <= 101 (supplement, sec:angular).
usage: python3 v6_light_action.py g N"""
import sys
from fractions import Fraction as Fr
import v3_lp_hp as V
from common import ladder_exact, W, allowed
g = int(sys.argv[1]); N = int(sys.argv[2])
o = V.run(g, N)
mus, lam = o['mus'], Fr(o['c2']) / 3
print(f"g={g} N={N}: |c|^2 = {o['c2']:.6e}")
worst = {}
for tw in [Fr(2), Fr(201, 100), Fr(3), Fr(4), Fr(6), Fr(g, 2), Fr(g) - Fr(1, 10)]:
    best = (0, None)
    for J in list(range(0, 4 * N + 8)) + [200, 500]:
        for b in (0, 4):
            if not allowed(b, J):
                continue
            D = J + tw
            r = ladder_exact(b, J, D, N)
            f = sum(m * x for m, x in zip(mus, r))
            val = (1 if b == 4 else -1) * f - lam * W(b, D)       # excess per unit |a_0|^2
            if val > best[0]:
                best = (val, (b, J))
    print(f"   twist {float(tw):7.3f}: max_(b,J) (sigma f - lam W) per unit |a_0|^2 = {float(best[0]):.3e} at {best[1]}")
