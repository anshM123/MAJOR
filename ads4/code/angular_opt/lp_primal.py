"""Primal LP: maximize |c|^2 subject to the endpoint rule and the ANEC harmonic-balance rules R_k (k=0..n).
States: (b, J, Delta) with twist >= g; b=0 for J in {0,2,4,6,...}; b=4 for J>=4 (all J>=4).
Variables: y = rho*W (share of the endpoint budget 3)."""
import numpy as np, sys
from scipy.optimize import linprog

def ladder(b, J, D, n):
    C = J*(J+1.0); T = (D-1.0)*(D-2.0)
    r = np.zeros(n+1); r[0] = 1.0
    for k in range(n):
        s = b + k
        if k+1 > J - b: break
        r[k+1] = r[k]*((k+4.0)/(k+1.0))**2*(C - s*(s+1))/(T - s*(s+1))
    return r

def W(b, D):
    T = (D-1.0)*(D-2.0)
    return (T - b*(b-1))/3.0 + 1.5

def states(g, Jmax, nt=40):
    out = []
    ts = np.concatenate([[0.0], np.geomspace(0.01*g, 50*g, nt)])
    for J in range(0, Jmax+1):
        for b in (0, 4):
            if b == 0 and J % 2 == 1: continue
            if b == 4 and J < 4: continue
            for t in ts:
                out.append((b, J, J + g + t))
    return out

def solve(g, n, Jmax, nt=40, verbose=False):
    S = states(g, Jmax, nt)
    cols = []
    for (b, J, D) in S:
        r = ladder(b, J, D, n)/W(b, D)
        cols.append(r if b == 4 else -r)
    A = np.array(cols).T   # (n+1) x N ; row k: m4(k)-m0(k)
    # objective: maximize row 2; constraints rows k!=2 equal 0; sum y <= 3
    rows = [k for k in range(n+1) if k != 2]
    sc = np.abs(A).max(axis=1); sc[sc == 0] = 1
    Aeq = A[rows]/sc[rows][:, None]
    beq = np.zeros(len(rows))
    res = linprog(-A[2], A_ub=np.ones((1, A.shape[1])), b_ub=[3.0], A_eq=Aeq, b_eq=beq,
                  bounds=(0, None), method='highs')
    c2 = -res.fun
    if verbose:
        y = res.x; idx = np.argsort(-y)[:12]
        for i in idx:
            if y[i] > 1e-12: print('   ', S[i], y[i])
    return c2, res

if __name__ == '__main__':
    for g in [10, 20, 40, 80]:
        for n in [2, 4, 8, 16, 32, 64]:
            c2, res = solve(g, n, Jmax=int(6*g)+10)
            tau = 4*np.sqrt(c2)
            print(f"g={g} n={n} |c|^2={c2:.4e} tau={tau:.4e} tau*g={tau*g:.3f} tau*g^3={tau*g**3:.1f} status={res.status}", flush=True)
