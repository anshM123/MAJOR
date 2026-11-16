"""
Linear programme: maximise |c|^2 subject to
   (S_k)  delta_{k,2}|c|^2 + sum_i p_i W0_k(Delta_i,J_i) = sum_j q_j W4_k(Delta_j,J_j) =: a_{2+k},   k = 0..kmax
   p_i, q_j >= 0;  EEC(theta) = 1 + 2 sum_{n>=2} a_n cos(n theta) >= 0 on a theta grid,
with the '++' data p_i on even spins J >= 0 and the '+-' data q_j on spins J >= 4, where every operator of spin > 2
has Delta >= max(gap, J + 1 + eps).  (Necessary conditions only: commutativity + positivity in the state |T,2>.)
A feasible point with spins <= Jmax is an exact solution of the infinite system (all sums are finite).
usage: lp_balance.py gap Jmax [xmax]     (xmax: optional upper bound on J/Delta for the spin >= 4 operators)
"""
import sys
import numpy as np, mpmath as mp
from scipy.optimize import linprog
import closed_forms as cf
mp.mp.dps = 30

def profiles(gap, Jmax, xmax=None, eps=1e-3):
    offs = [0, 0.02, 0.1, 0.3, 0.6, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 256, 512, 2048]
    kmax = Jmax
    pp, pm = [], []
    def add(lst, W, J, top):
        base = max(gap, J + 1 + eps) if J > 2 else (J + 1 + eps if J > 0 else 0.5 + eps)
        if xmax is not None and J > 2:
            base = max(base, J/xmax)
        for o in offs:
            Delta = mp.mpf(base + o)
            v = np.array([float(W(k, Delta, J)) for k in range(kmax + 1)])
            if not np.all(np.isfinite(v)) or v.max() <= 0:
                continue
            sc = v.max()
            lst.append((J, float(Delta), v/sc, sc))
            if J == 0:
                break
    for J in range(0, Jmax + 1, 2):
        add(pp, cf.W0, J, J)
    for J in range(4, Jmax + 5):
        add(pm, cf.W4, J, J - 4)
    return pp, pm, kmax

def solve(gap, Jmax, xmax=None, ntheta=721, verbose=True):
    pp, pm, kmax = profiles(gap, Jmax, xmax)
    npp, npm = len(pp), len(pm)
    nv = 1 + npp + npm                      # variables: c2, p (rescaled), q (rescaled)
    Aeq = np.zeros((kmax + 1, nv)); beq = np.zeros(kmax + 1)
    for k in range(kmax + 1):
        if k == 2:
            Aeq[k, 0] = 1
        for i, (J, D, v, sc) in enumerate(pp):
            Aeq[k, 1 + i] = v[k]
        for j, (J, D, v, sc) in enumerate(pm):
            Aeq[k, 1 + npp + j] = -v[k]
    # EEC >= 0:  -2 sum_k a_{2+k} cos((2+k) theta) <= 1, a_{2+k} = sum_j q_j W4_k
    th = np.linspace(0, np.pi, ntheta)
    Aub = np.zeros((ntheta, nv)); bub = np.ones(ntheta)
    for j, (J, D, v, sc) in enumerate(pm):
        Aub[:, 1 + npp + j] = -2*sum(v[k]*np.cos((2 + k)*th) for k in range(kmax + 1))
    cost = np.zeros(nv); cost[0] = -1
    res = linprog(cost, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=[(0, None)]*nv, method='highs')
    if verbose:
        print(f"gap={gap}, Jmax={Jmax}, xmax={xmax}: status={res.status} max |c|^2 = {-res.fun if res.status == 0 else None}")
        if res.status == 0:
            x = res.x
            a = np.array([sum(x[1 + npp + j]*pm[j][2][k] for j in range(npm)) for k in range(kmax + 1)])
            print("   a_{2+k}, k=0..9:", np.round(a[:10], 4))
            sup = [(pp[i][0], round(pp[i][1], 3), x[1 + i]/pp[i][3]) for i in range(npp) if x[1 + i] > 1e-9]
            sum_ = [(pm[j][0], round(pm[j][1], 3), x[1 + npp + j]/pm[j][3]) for j in range(npm) if x[1 + npp + j] > 1e-9]
            print("   '++' operators (J, Delta, |a^{++}|^2):", [(J, D, float(f"{w:.3g}")) for J, D, w in sup][:14], '...' if len(sup) > 14 else '')
            print("   '+-' operators (J, Delta, |a^{+-}|^2):", [(J, D, float(f"{w:.3g}")) for J, D, w in sum_][:14], '...' if len(sum_) > 14 else '')
            resid = np.abs(Aeq @ x - beq).max()
            print(f"   max residual of the sum rules: {resid:.2e};  min EEC on grid: {(1 - Aub @ x).min():.3e}")
    return res, pp, pm

if __name__ == '__main__':
    gap = float(sys.argv[1]); Jmax = int(sys.argv[2])
    xmax = float(sys.argv[3]) if len(sys.argv) > 3 else None
    solve(gap, Jmax, xmax)
