"""Triage: quantum Mrs. Gerber lower bound (Hirche-Reeb conjecture).
X_i uniform bits with quantum side information rho_x^{(i)} on B_i; Z = X1 xor X2; s_i = H(X_i|B_i) (bits).
Conjecture: H(Z|B1B2) >= F(s1,s2) if s1+s2<=1, else s1+s2-1+F(1-s1,1-s2),  F(s,t)=h(h^{-1}(s)*h^{-1}(t))."""
import numpy as np, sys
from scipy.optimize import minimize, brentq
def h2(p): 
    p = min(max(p, 0.0), 1.0)
    return 0.0 if p in (0.0, 1.0) else -p*np.log2(p)-(1-p)*np.log2(1-p)
def h2inv(s):
    s = min(max(s, 0.0), 1.0)
    if s <= 0: return 0.0
    if s >= 1: return 0.5
    return brentq(lambda p: h2(p)-s, 0, 0.5, xtol=1e-15)
def F(s, t):
    a, b = h2inv(s), h2inv(t); return h2(a*(1-b)+(1-a)*b)
def bound(s1, s2): return F(s1, s2) if s1+s2 <= 1 else s1+s2-1+F(1-s1, 1-s2)
def S(M):
    w = np.linalg.eigvalsh((M+M.conj().T)/2); w = w[w > 1e-300]; return float(-(w*np.log2(w)).sum())
def st(x, d, r):
    L = (x[:d*r]+1j*x[d*r:]).reshape(d, r); M = L@L.conj().T; return M/np.trace(M).real
def quantities(x, d1, d2, r):
    k1 = 2*d1*r; k2 = 2*d2*r
    a0, a1 = st(x[:k1], d1, r), st(x[k1:2*k1], d1, r); b0, b1 = st(x[2*k1:2*k1+k2], d2, r), st(x[2*k1+k2:], d2, r)
    s1 = 1 + 0.5*(S(a0)+S(a1)) - S((a0+a1)/2); s2 = 1 + 0.5*(S(b0)+S(b1)) - S((b0+b1)/2)
    g0 = 0.5*(np.kron(a0, b0)+np.kron(a1, b1)); g1 = 0.5*(np.kron(a0, b1)+np.kron(a1, b0))
    HZ = 1 + 0.5*(S(g0)+S(g1)) - S((g0+g1)/2)
    return HZ, s1, s2
def gap(x, d1, d2, r):
    HZ, s1, s2 = quantities(x, d1, d2, r); return HZ - bound(s1, s2)
if __name__ == "__main__":
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    for (d1, d2, r) in [(2, 2, 1), (2, 2, 2), (3, 3, 1), (2, 3, 2), (3, 3, 3), (4, 4, 1), (4, 4, 2)]:
        best = (np.inf, None)
        for t in range(40):
            x0 = rng.normal(size=4*d1*r+4*d2*r)
            res = minimize(gap, x0, args=(d1, d2, r), method="L-BFGS-B", options=dict(maxiter=300))
            if res.fun < best[0]: best = (res.fun, res.x)
        HZ, s1, s2 = quantities(best[1], d1, d2, r)
        print(f"dB=({d1},{d2}) rank={r}: min [H(Z|B1B2) - bound] = {best[0]:+.3e}  at s1={s1:.4f} s2={s2:.4f} HZ={HZ:.6f}" + ("   <<< VIOLATION" if best[0] < -1e-7 else ""), flush=True)
