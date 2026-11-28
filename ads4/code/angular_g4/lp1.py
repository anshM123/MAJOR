"""LP1: primal  max |c|^2  s.t. endpoint (=3), R_k (k<=K), M_n (2<=n<=N) over a grid of heavy states (twist>=g).
Valid as an upper bound only through its dual (finite rule set = relaxation); primal on a finite grid is a lower
estimate of the true LP value of that rule set."""
import sys, numpy as np
from scipy.optimize import linprog
from core import ladder, hvec, W, allowed, GAMMA


def grid(g, Jmax_lin, Jmax, nJlog=80, nt=60, tmax=1e4):
    Js = sorted(set(list(range(0, Jmax_lin + 1)) + [int(x) for x in np.geomspace(Jmax_lin + 1, Jmax, nJlog)]))
    ts = np.concatenate([[0.0], np.geomspace(1e-3, tmax, nt) * g])
    return [(b, J, J + g + t) for J in Js for b in (0, 4) if allowed(b, J) for t in ts]


def columns(S, K, N):
    cols = []
    for (b, J, D) in S:
        r = ladder(b, J, D, max(K, N + 1))
        h = hvec(b, J, D, N)
        w = W(b, D)
        sig = 1.0 if b == 4 else -1.0
        cols.append(np.concatenate([[w], sig * r[:K + 1], h[2:N + 1]]) / w)   # variable y = rho*W
    return np.array(cols).T


def solve(g, K, N, S=None, use_M=True, use_end=True, verbose=0):
    if S is None:
        S = grid(g, 12 * g, 2000 * g)
    A = columns(S, K, N)
    nv = A.shape[1]
    # extra variable c2 (last)
    rows = []; rhs = []
    if use_end:
        rows.append(np.concatenate([A[0], [0.0]])); rhs.append(3.0)
    for k in range(K + 1):
        rows.append(np.concatenate([A[1 + k], [-1.0 if k == 2 else 0.0]])); rhs.append(0.0)
    if use_M:
        for n in range(2, N + 1):
            rows.append(np.concatenate([A[K + 2 + n - 2], [GAMMA.get(n, 0.0)]])); rhs.append(1.0)
    Aeq = np.array(rows); beq = np.array(rhs)
    sc = np.abs(Aeq).max(axis=1); sc[sc == 0] = 1
    c = np.zeros(nv + 1); c[-1] = -1.0
    res = linprog(c, A_eq=Aeq / sc[:, None], b_eq=beq / sc, bounds=[(0, None)] * (nv + 1), method='highs')
    if res.status != 0:
        return None, res, S
    if verbose:
        y = res.x[:-1]; idx = np.argsort(-y)[:verbose]
        for i in idx:
            if y[i] > 1e-14: print('    ', S[i][0], S[i][1], round(S[i][2] - S[i][1] - g, 4), '%.3e' % y[i])
    return res.x[-1], res, S


if __name__ == '__main__':
    g = int(sys.argv[1]); K = int(sys.argv[2]); N = int(sys.argv[3])
    S = grid(g, 12 * g, 2000 * g, nJlog=60, nt=40)
    print('states', len(S))
    for useM in (False, True):
        c2, res, _ = solve(g, K, N, S, use_M=useM, verbose=12)
        if c2 is None:
            print('useM', useM, 'status', res.status, res.message); continue
        tau = 4 * np.sqrt(max(c2, 0))
        print(f'g={g} K={K} N={N} useM={useM}: |c|^2={c2:.5e} tau={tau:.5e} tau g^3={tau*g**3:.2f} tau g^4={tau*g**4:.1f}')
