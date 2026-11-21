"""
po_pwtest (PARITY_ODD_CLOSE.md sect. 1): POINT-WISE test of the OO <-> ++ identity on the 3x3 [even(2), odd(1)] integrand of the
EE functional (w0_ee_common.make_build_odd_EE on p8's point-wise copy of p3_odd_action) and of C^_3 (p3's own o_H q, kind C3H).
At random points z = z* + L xi of the 7-dim integration space (n0 = +1 and n0 = -1), the prefactor matrix pref[a, c](z) (a: 12-side
structure, c: 34-side structure; numerators x structure denominators, before the shadow inverse) is tested for
  (R1) row relation:   pref[O, c](z) = sum_a w_a pref[a, c](z)   for all c, z, with CONSTANT w (12 side)
  (R2) column relation: pref[a, O](z) = sum_c v_c pref[a, c](z)  (34 side)
  (R3) OO = (w^T pref_EE v)  pointwise.
Read-only imports. Usage: python3 po_pwtest.py J Delta nu [kind] [scale]
"""
import os, sys
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6'); os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import l1_tps as T
import l1_crystal as Cr
import p8_odd_pw as P8M
import n1_pw as PW
import w0_ee_common as W0

P8 = P8M.P8
BEE = W0.make_build_odd_EE(P8)


def saddle_L(nu, Delta, n0=1.0):
    zs, a = Cr.saddle(nu, Delta)
    zv = [T.var(i, zs[i]) for i in range(T.NV)]
    lamF = Cr.build(zv, nu, Delta, 0, n0=n0)[0]
    _, HF = T.grad_hess(lamF)
    L = np.linalg.inv(T.csym_cholesky(-HF)).T
    for _ in range(2):
        lamF2 = Cr.build(T.linsub(zs, L), nu, Delta, 0, n0=n0)[0]
        _, H2 = T.grad_hess(lamF2)
        L = L @ np.linalg.inv(T.csym_cholesky(-H2)).T
    return np.asarray(zs, complex), L


def prefs(J, Delta, nu, kind='EE', n0=1.0, npts=40, scale=0.5, seed=1):
    zs, L = saddle_L(nu, Delta, n0)
    rng = np.random.default_rng(seed)
    xi = scale*rng.standard_normal((npts, 7))
    zl = PW.points(zs, L, xi)
    if kind == 'EE':
        out = BEE(zl, nu, Delta, J, 'EE', n0)
    else:
        out = P8.build(zl, nu, Delta, J, kind, n0)
    return out[2]                                  # (na, nc, N)


def lsq_rel(A, b):
    """b (n) ~ A (n x k) w; returns w, relative residual."""
    w, *_ = np.linalg.lstsq(A, b, rcond=None)
    return w, np.linalg.norm(A @ w - b)/np.linalg.norm(b)


if __name__ == '__main__':
    J, Dl, nu = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    kind = sys.argv[4] if len(sys.argv) > 4 else 'EE'
    sc = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
    for n0 in (1.0, -1.0):
        P = prefs(J, Dl, nu, kind, n0, scale=sc)
        na = P.shape[0]
        ev = list(range(na - 1)); O = na - 1
        P = P/np.abs(P[ev][:, ev]).max(axis=(0, 1))[None, None, :]     # per-point normalisation (common factor)
        A = np.concatenate([P[ev, c, :].T for c in range(na)], axis=0)  # rows: (c, point), cols: a in ev
        b = np.concatenate([P[O, c, :] for c in range(na)])
        w, r1 = lsq_rel(A, b)
        A2 = np.concatenate([P[a, ev, :].T for a in range(na)], axis=0)
        b2 = np.concatenate([P[a, O, :] for a in range(na)])
        v, r2 = lsq_rel(A2, b2)
        oo_pred = np.einsum('a,c,acn->n', w, v, P[np.ix_(ev, ev)])
        r3 = np.linalg.norm(oo_pred - P[O, O])/np.linalg.norm(P[O, O])
        print(f"J={J} D={Dl} nu={nu} kind={kind} n0={n0:+.0f}: R1 w={np.round(w, 6)} resid {r1:.2e} | R2 v={np.round(v, 6)} "
              f"resid {r2:.2e} | R3 resid {r3:.2e} | |EO|/|EE| {np.abs(P[ev, O]).mean()/np.abs(P[np.ix_(ev, ev)]).mean():.3f}",
              flush=True)
