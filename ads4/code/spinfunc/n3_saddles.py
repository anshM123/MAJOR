"""
n3_saddles (doc 25): the critical points of the scalar flow function Phi_nu (n3_flow / n1_direct.phi_w) in the soft block at small nu.

Finding [NE]: the crystal saddle S (symmetric under sigma: (w1, w2, w3, w4) -> conj(w3, w4, w1, w2), w0 = w5 = w6 = 0) undergoes a
pitchfork at nu* (det of the soft 4x4 block -> 0; nu* = 11.28 at Delta = 100, J = 4); below nu* a sigma-conjugate pair A+-, A- = sigma(A+)
branches off.  find(nu, Delta, J, which) continues each branch from nu = 20 downwards (S: Newton + symmetrisation after each step).
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import n1_direct as D
import l1_tps as T
import l1_crystal as Cr
import n3_flow as F


def sigma(w):
    return np.array([0.0, np.conj(w[3]), np.conj(w[4]), np.conj(w[1]), np.conj(w[2]), 0.0, 0.0], complex)


def _gh(w, c, nu, Delta, J):
    wv = [T.var(i, w[i]) for i in range(T.NV)]
    P = D.phi_w(wv, c, nu, Delta, J)
    g, H = T.grad_hess(P)
    return g, H, P[0]


def newton(w, nu, Delta, J, sym=False, it=60, cap=0.05, tol=1e-12):
    c, _ = Cr.saddle(nu, Delta)
    for k in range(it):
        g, H, _ = _gh(w, c, nu, Delta, J)
        dw = np.linalg.solve(H, g)
        s = 1.0 if np.abs(dw).max() <= cap else cap/np.abs(dw).max()
        w = w - s*dw
        if sym:
            w = 0.5*(w + sigma(w))
        if np.abs(dw).max() < tol:
            break
    g, H, P0 = _gh(w, c, nu, Delta, J)
    return w, g, H, P0, c


_C = {}


def find(nu, Delta, J, which='S', nu0=20.0, step=0.25):
    """returns (c, wc, L, info) for the branch which in {'S', 'A+', 'A-'} at nu (continuation from nu0)."""
    key = (round(nu, 9), float(Delta), int(J), which)
    if key in _C:
        return _C[key]
    base = 'S' if which == 'S' else 'A'
    w = np.zeros(7, complex)
    nus = list(np.arange(nu0, nu, -step)) + [nu]
    for v in nus:
        w, g, H, P0, c = newton(w, float(v), Delta, J, sym=(base == 'S'))
    if which == 'A-':
        w, g, H, P0, c = newton(sigma(w), nu, Delta, J)
    elif which == 'A+':
        pass
    asym = abs(w[1] - np.conj(w[3])) + abs(w[2] - np.conj(w[4]))
    L = np.linalg.inv(T.csym_cholesky(-H)).T
    info = dict(grad=float(np.abs(g).max()), asym=float(asym), Phi0=P0, detsoft=complex(np.linalg.det(-H[1:5, 1:5])),
                ev=np.sort_complex(np.linalg.eigvals(-H)))
    _C[key] = (c, w, L, info)
    return _C[key]


if __name__ == '__main__':
    nu, Dl, J = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
    for wh in ('S', 'A+', 'A-'):
        c, w, L, info = find(nu, Dl, J, wh)
        print(wh, np.round(w, 4), {k: (np.round(v, 4) if not isinstance(v, np.ndarray) else np.round(v, 2)) for k, v in info.items()}, flush=True)
