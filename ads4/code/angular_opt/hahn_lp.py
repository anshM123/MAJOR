"""LP for the Lemma-M certificate: q polynomial of degree <= n on the grid {0..M}, q(0)=1, -theta<=q<=1 on the grid;
maximise c2 = C(M,2)(q(2)-2q(1)+q(0)).  Basis: discrete Chebyshev (Gram) polynomials, orthonormal on {0..M}."""
import numpy as np, sys
from scipy.optimize import linprog
def gram_basis(M, n):
    z = np.arange(M+1, dtype=float)
    P = np.zeros((n+1, M+1)); P[0] = 1.0/np.sqrt(M+1)
    if n >= 1:
        v = z - z.mean(); P[1] = v/np.linalg.norm(v)
    for j in range(1, n):
        # Stieltjes / Lanczos three-term recurrence with full reorthogonalisation
        w = (z - z.mean())*P[j]
        for i in range(j+1): w -= np.dot(w, P[i])*P[i]
        for i in range(j+1): w -= np.dot(w, P[i])*P[i]
        P[j+1] = w/np.linalg.norm(w)
    return P
def solve(M, n, theta):
    G = gram_basis(M, n)            # (n+1) x (M+1)
    obj = -(G[:, 2] - 2*G[:, 1])
    Aub = np.vstack([G.T, -G.T]); bub = np.concatenate([np.ones(M+1), theta*np.ones(M+1)])
    Aeq = G[:, 0:1].T; beq = [1.0]
    res = linprog(obj, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=[(None, None)]*(n+1), method='highs')
    a = res.x; q = a @ G
    c2 = M*(M-1)/2*(q[2]-2*q[1]+q[0])
    return c2, q, a, res
if __name__ == '__main__':
    g = int(sys.argv[1]); M = g-2
    T = lambda D: (D-1)*(D-2)
    WB = (T(g+4)-12)/3+1.5; W0 = T(g+2)/3+1.5; theta = W0/WB
    for n in map(int, sys.argv[2:]):
        c2, q, a, res = solve(M, n, theta)
        L = WB*c2/100; tau = 4*np.sqrt(3/L)
        print(f"g={g} M={M} n={n} status={res.status} c2/M^2={c2/M**2:.4f} L={L:.4e} tau={tau:.4e} tau*g^2={tau*g*g:.2f} q1={q[1]:.4f} q2={q[2]:.4f}", flush=True)
