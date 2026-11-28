"""Dual LP: maximize L=-lambda_2 s.t. for every heavy state
   b=4 (J>=4):        W4 + sum_k lam_k r_k^(4) >= 0
   b=0 (J even >=0):  W0 - sum_k lam_k r_k^(0) >= 0
Then |c|^2 <= (3+D)/L, tau^2 = 16|c|^2."""
import numpy as np
from scipy.optimize import linprog
from lp_primal import ladder, W, states

def build(S, n):
    rows = []
    for (b, J, D) in S:
        r = ladder(b, J, D, n)/W(b, D)
        rows.append(r if b == 4 else -r)   # constraint: 1 + row.lam >= 0
    return np.array(rows)

def solve_dual(g, n, Jmax, nt=40):
    S = states(g, Jmax, nt)
    R = build(S, n)
    sc = np.abs(R).max(axis=0); sc[sc == 0] = 1
    Rs = R/sc
    # variables mu_k = lam_k*sc_k ; objective maximize -lam_2 = -mu_2/sc_2
    c = np.zeros(n+1); c[2] = 1.0/sc[2]
    res = linprog(c, A_ub=-Rs, b_ub=np.ones(len(S)), bounds=[(None, None)]*(n+1), method='highs')
    lam = res.x/sc
    return -lam[2], lam, res

def check(g, lam, Jmax, nt=400):
    n = len(lam)-1
    worst = (1e9, None)
    ts = np.concatenate([[0.0], np.geomspace(1e-4*g, 1e4*g, nt)])
    for J in range(0, Jmax+1):
        for b in (0, 4):
            if b == 0 and J % 2: continue
            if b == 4 and J < 4: continue
            for t in ts:
                D = J+g+t
                r = ladder(b, J, D, n)
                v = 1 + (1 if b == 4 else -1)*np.dot(lam, r)/W(b, D)
                if v < worst[0]: worst = (v, (b, J, D))
    return worst

if __name__ == '__main__':
    for g in [10, 20]:
        L, lam, res = solve_dual(g, int(1.2*g), 6*g+10)
        print(g, 'L=', L, 'L/g^6*7200=', L/g**6*7200, 'tau=', 4*np.sqrt(3/L), res.status)
        print('  lam:', lam[:8])
        print('  check up to J=50g:', check(g, lam, 50*g, nt=200))
