"""Dual LP with LIGHT-state robustness: maximise L=-lambda_2 subject to
   heavy (twist>=g): W4 + Lam4 >= 0 ; W0 - Lam0 >= 0
   light (2 < twist < g): |Lam_b| <= kappa * W_b
so that the light-state term in the bound is <= 16 kappa E'_light / L, with E'_light the light share of the endpoint rule."""
import numpy as np, sys
from scipy.optimize import linprog
from lp_primal import ladder, W

def rows(g, n, kappa):
    A = []; bnd = []
    Js = sorted(set(list(range(0, 6*g)) + [int(x) for x in np.geomspace(6*g, 3000*g, 120)]))
    ts_h = np.concatenate([[0.0], np.geomspace(0.003*g, 300*g, 30)])
    for J in Js:
        for b in (0, 4):
            if b == 0 and J % 2: continue
            if b == 4 and J < 4: continue
            for t in ts_h:
                D = J + g + t; r = ladder(b, J, D, n)/W(b, D)
                A.append(-r if b == 4 else r); bnd.append(1.0)      # -(lam.r)/W <= 1  (b4) ; (lam.r)/W <= 1 (b0)
            if kappa is None: continue
            for h in np.concatenate([[2.0 + 1e-9, 2.05, 2.3], np.linspace(2.6, g-1e-6, 25)]):
                D = J + h; r = ladder(b, J, D, n)/W(b, D)
                A.append(r); bnd.append(kappa); A.append(-r); bnd.append(kappa)
    return np.array(A), np.array(bnd)

def solve(g, n, kappa):
    A, bnd = rows(g, n, kappa)
    sc = np.abs(A).max(axis=0); sc[sc == 0] = 1
    c = np.zeros(n+1); c[2] = 1/sc[2]
    res = linprog(c, A_ub=A/sc, b_ub=bnd, bounds=[(None, None)]*(n+1), method='highs')
    if res.status != 0: return None, res
    lam = res.x/sc
    return -lam[2], lam

if __name__ == '__main__':
    g = int(sys.argv[1]); n = int(sys.argv[2])
    for kappa in [None, 1e8, 1e6, 1e4, float(g)**3, float(g)**2, 10*g, 1.0]:
        L, lam = solve(g, n, kappa)
        if L is None: print(g, n, kappa, 'fail', lam.message); continue
        tau = 4*np.sqrt(3/L)
        print(f"g={g} n={n} kappa={kappa}: L={L:.4e} tau={tau:.4e} tau*g={tau*g:.3f} tau*g^2={tau*g*g:.1f}", flush=True)
