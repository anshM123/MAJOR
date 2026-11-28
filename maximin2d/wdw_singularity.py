"""Does the Wheeler-DeWitt set W(C) of C = {t_L = t_R = 0} reach the singularity in thin-shell
Vaidya-BTZ (shell from the right boundary at v0 < 0, mass M1 -> M2)?  Ingoing EF (right chart):
ds^2 = -(r^2 - M(v))dv^2 + 2 dv dr + r^2 dphi^2.  Top vertex of W(C) = meeting point of the right ingoing
ray v = 0 and the left ray from t_L = 0, which inside the future interior obeys dr/dv = (r^2 - M(v))/2.
In pure BTZ(M1) the left ray hits r = 0 exactly at v = 0 (square Penrose diagram).  Here: integrate the
pure-M1 trajectory up to v0, then M2.  If r reaches 0 at v* < 0, W(C) meets the singularity (K fails).
usage: python3 -I wdw_singularity.py"""
import numpy as np
from scipy.integrate import solve_ivp
M1 = 1.0
m = np.sqrt(M1)
# pure BTZ(M1) left ray reaching r=0 at v=0: dr/dv=(r^2-1)/2 -> r(v) = tanh(-v/2) (check: r(0)=0, r->1 as v->-inf)
r_M1 = lambda v: np.tanh(-m * v / 2) * m
for v0 in (-0.5, -1.0, -2.0):
    for M2 in (1.5, 1.1, 0.9, 0.5):
        def rhs(v, r): return [(r[0] ** 2 - M2) / 2]
        def hit(v, r): return r[0]
        hit.terminal = True
        sol = solve_ivp(rhs, (v0, 5.0), [r_M1(v0)], events=hit, rtol=1e-12, atol=1e-14)
        vs = sol.t_events[0][0] if sol.t_events[0].size else np.inf
        tag = "NEC" if M2 > M1 else "NEC-viol"
        print(f"v0={v0:+.1f} M2={M2:.1f} ({tag:8s}): left ray reaches r=0 at v*={vs:+.4f} ->",
              "W(C) meets singularity (K fails)" if vs < 0 else "top vertex in the interior (K ok here)")
