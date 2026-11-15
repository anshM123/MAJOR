"""Triage: two-copy multiplicativity of the maximal output p-norm for polarized Werner-Holevo channels,
Phi_x = x id + (1-x) W_d,  W_d(X) = (tr(X) 1 - X^T)/(d-1),  1<p<2, d>=3."""
import numpy as np, sys
from scipy.optimize import minimize
def chan(rho, d, x): return x*rho + (1-x)*(np.trace(rho)*np.eye(d) - rho.T)/(d-1)
def pnorm(M, p):
    w = np.linalg.eigvalsh((M+M.conj().T)/2); return (np.abs(w)**p).sum()**(1/p)
def one(v, d, x, p):
    psi = v[:d] + 1j*v[d:]; psi /= np.linalg.norm(psi); return -pnorm(chan(np.outer(psi, psi.conj()), d, x), p)
def two_out(rho, d, x):
    # apply Phi (x) Phi to a d^2 x d^2 matrix
    T = rho.reshape(d, d, d, d)                      # (a b ; a' b')
    def ap1(T):  # on first factor
        tr = np.einsum('abac->bc', T); out = x*T + (1-x)*(np.einsum('ac,bd->abcd', np.eye(d), tr) - T.transpose(2, 1, 0, 3))/(d-1); return out
    def ap2(T):
        tr = np.einsum('abcb->ac', T); out = x*T + (1-x)*(np.einsum('bd,ac->abcd', np.eye(d), tr) - T.transpose(0, 3, 2, 1))/(d-1); return out
    return ap2(ap1(T)).reshape(d*d, d*d)
def two(v, d, x, p):
    psi = v[:d*d] + 1j*v[d*d:]; psi /= np.linalg.norm(psi); return -pnorm(two_out(np.outer(psi, psi.conj()), d, x), p)
if __name__ == "__main__":
    rng = np.random.default_rng(1)
    # self-test of two_out on product input
    d = 3; a = rng.normal(size=d)+1j*rng.normal(size=d); b = rng.normal(size=d)+1j*rng.normal(size=d)
    ra = np.outer(a, a.conj()); rb = np.outer(b, b.conj())
    assert np.allclose(two_out(np.kron(ra, rb), d, .3), np.kron(chan(ra, d, .3), chan(rb, d, .3)))
    for d in (3, 4):
        for p in (1.05, 1.2, 1.5, 1.8, 1.95):
            for x in (0.05, 0.15, 0.3, 0.5, 0.7, 0.9):
                n1 = max(-minimize(one, rng.normal(size=2*d), args=(d, x, p), method="BFGS").fun for _ in range(8))
                n2 = max(-minimize(two, rng.normal(size=2*d*d), args=(d, x, p), method="L-BFGS-B", options=dict(maxiter=400)).fun for _ in range(12))
                # also try the maximally entangled input explicitly
                me = np.eye(d).reshape(-1)/np.sqrt(d); n2 = max(n2, -two(np.concatenate([me, 0*me]), d, x, p))
                print(f"d={d} p={p} x={x}: nu_p^2={n1**2:.8f}  two-copy max={n2:.8f}  ratio={n2/n1**2:.8f}" + ("   <<< VIOLATION" if n2 > n1**2*(1+1e-7) else ""), flush=True)
