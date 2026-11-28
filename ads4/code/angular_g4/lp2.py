"""Dual LP with endpoint (multiplier nu>=0), R_k (lambda_k, k<=K), per-harmonic M_n (mu_n, 2<=n<=N).
Heavy:  F = nu W + sum mu_n h_n + sigma sum lambda_k r_k >= 0.
Light (2<twist<g), optional: mode 'theta': F >= -Theta W  (light coef 16 Theta/L);  mode 'nu': F >= -(kappa-1) nu W.
Normalisation 3 nu + sum mu_n = 3.  Maximise L = sum mu_n gamma_n - lambda_2.
Then  |c|^2 <= (3 + (kappa-1) nu E'_light)/L,   tau^2 = 16|c|^2.   (D = 0.)"""
import numpy as np, sys, time
from scipy.optimize import linprog
from vec import rmat, hmat, Wv, states

GAM = {3: -0.5, 4: 0.5, 5: 1.5}


def build(B, J, D, K, N):
    KK = max(K, N + 1)
    r = rmat(B, J, D, KK)
    h = hmat(B, J, D, N, r) if N >= 2 else np.zeros((len(J), 1))
    w = Wv(B, D)
    sig = np.where(B == 4, 1.0, -1.0)
    # columns: nu, mu_2..mu_N, lam_0..lam_K ; each row divided by w
    cols = [np.ones(len(J))]
    for n in range(2, N + 1):
        cols.append(h[:, n] / w)
    for k in range(K + 1):
        cols.append(sig * r[:, k] / w)
    return np.array(cols).T


def solve(g, K, N, kappa=None, mode='theta', Jlin=None, Jmax=None, nJlog=80, nt=50, nlight=40, ret=False):
    Jlin = Jlin or 12 * g; Jmax = Jmax or 3000 * g
    ts = np.concatenate([[0.0], np.geomspace(1e-3, 1e4, nt) * g])
    lh = None
    if kappa is not None:
        lh = np.concatenate([[2.0 + 1e-9, 2.01, 2.05, 2.2, 2.5, 3, 4, 5, 6, 8], np.geomspace(10, g - 1e-6, nlight)])
    B, J, D, Lt = states(g, Jlin, Jmax, nJlog, ts, lh)
    A = build(B, J, D, K, N)
    nM = max(N - 1, 0)
    # inequality A x >= 0 (heavy) ; light: A x + (kappa-1) nu >= 0
    Aub = -A.copy()
    if kappa is not None:
        if mode == 'nu':
            Aub[Lt == 1, 0] -= (kappa - 1.0)
    bub = np.zeros(len(J))
    if kappa is not None and mode == 'theta':
        bub[Lt == 1] = kappa          # F_light/W >= -Theta
    sc = np.abs(Aub).max(axis=0); sc[sc == 0] = 1
    nv = A.shape[1]
    # objective: maximise L = sum mu_n gamma_n - lambda_2
    c = np.zeros(nv)
    for n in range(2, N + 1):
        c[1 + n - 2] = -GAM.get(n, 0.0)
    c[1 + nM + 2] = 1.0
    Aeq = np.zeros((1, nv)); Aeq[0, 0] = 3.0; Aeq[0, 1:1 + nM] = 1.0
    bounds = [(0, None)] + [(None, None)] * (nv - 1)
    res = linprog(c / sc, A_ub=Aub / sc, b_ub=bub/1.0, A_eq=Aeq / sc, b_eq=[3.0], bounds=bounds, method='highs')
    if res.status != 0:
        return None
    x = res.x / sc
    L = -res.fun
    out = dict(L=L, nu=x[0], mu=x[1:1 + nM], lam=x[1 + nM:], tau=4 * np.sqrt(3 / L) if L > 0 else np.inf)
    if ret:
        out.update(B=B, J=J, D=D, Lt=Lt, A=A)
    return out


if __name__ == '__main__':
    g = int(sys.argv[1]); kap = None if sys.argv[2] == 'inf' else float(eval(sys.argv[2]))
    for N in [int(x) for x in sys.argv[3].split(',')]:
        K = max(N, int(sys.argv[4]) if len(sys.argv) > 4 else 12)
        for useM in (False, True):
            t0 = time.time()
            o = solve(g, K, N if useM else 0, kap)
            if o is None: print(g, N, useM, 'fail', flush=True); continue
            lc = 16 * kap / o['L'] if kap else 0
            print(f"g={g} kappa={kap} K={K} N={N if useM else 0}: L={o['L']:.4e} L/g^2={o['L']/g**2:.3f} nu={o['nu']:.3f} tau_heavy={o['tau']:.4e} tau*g={o['tau']*g:.3f} tau*g^3={o['tau']*g**3:.1f} lightcoef={lc:.4g} ({time.time()-t0:.0f}s)", flush=True)
