"""Light-robust dual LP with sparser J sampling (for large g)."""
import numpy as np, sys
from scipy.optimize import linprog
from lp_primal import ladder, W
def rows(g, n, kappa, nJ=400, nth=30, nlight=30):
    A = []; bnd = []
    Js = sorted(set(list(range(0, 300)) + [int(x) for x in np.geomspace(300, 3000*g, nJ)]))
    ts_h = np.concatenate([[0.0], np.geomspace(0.001*g, 1000*g, nth)])
    hl = np.concatenate([[2.0+1e-9, 2.01, 2.05, 2.2, 2.5, 3, 4, 5, 6, 8], np.geomspace(10, g-1e-6, nlight)])
    for J in Js:
        for b in (0, 4):
            if b == 0 and J % 2: continue
            if b == 4 and J < 4: continue
            for t in ts_h:
                D = J + g + t; r = ladder(b, J, D, n)/W(b, D)
                A.append(-r if b == 4 else r); bnd.append(1.0)
            for h in hl:
                D = J + h; r = ladder(b, J, D, n)/W(b, D)
                A.append(r); bnd.append(kappa); A.append(-r); bnd.append(kappa)
    return np.array(A), np.array(bnd)
def solve(g, n, kappa, **kw):
    A, bnd = rows(g, n, kappa, **kw)
    sc = np.abs(A).max(axis=0); sc[sc == 0] = 1
    c = np.zeros(n+1); c[2] = 1/sc[2]
    res = linprog(c, A_ub=A/sc, b_ub=bnd, bounds=[(None, None)]*(n+1), method='highs')
    if res.status != 0: return None, res
    return -res.x[2]/sc[2], res.x/sc
if __name__ == '__main__':
    g = int(sys.argv[1]); n = int(sys.argv[2])
    K = 729/4800*(g-2)**2*(1-49/(2*g+9)**2)
    print('paper K =', K, ' paper tau =', 4*np.sqrt(3/K), ' paper light coef ~ 16*70 = 1120')
    for kappa in [float(x) for x in sys.argv[3:]]:
        L, lam = solve(g, n, kappa)
        if L is None: print(g, n, kappa, 'fail', flush=True); continue
        print(f"g={g} n={n} kappa={kappa:.3g}: L={L:.4e} L/K={L/K:.2f} tau_heavy={4*np.sqrt(3/L):.4e} lightcoef={16*(kappa-1)/L:.4g}", flush=True)
