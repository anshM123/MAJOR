"""
n3fix/nf_gamma.py (N3C3H_FIX.md §3.3, §4): Gamma (c2, kappa subtracted) and the (EDGE-R) costs of the J = 26 gap-mass state (26, 1.026)
with the LOCAL patched-n3 / n1 points of this directory injected into EDGE_R2's batch-3 row, by EDGE_R2's own chain (edger2/e2_gamma.py,
imported read-only and run in-process; every write redirected to n3fix/cache/gamma/).

Basis.  EDGE_R2 maps the external records into the chain's frame, E' = A_ext E_ext B_ext (e2_basis.fix).  Local records are in this
machine's <TTO> structure basis, E_ext = G^T E_loc H (BATCH5_FOLD §1.2), so E' = A' E_loc B' with constant 2 x 2 A', B' per Delta (state,
reference).  A', B' are fitted (alternating least squares, full integrals exp(logN) E, EE and C3H, mean + 8 replicas) on two ANCHORS:
batch-3 n1 nodes (x = 0.9198, 0.8572) recomputed locally with the same n1 code path (cache/anchor).  Every injected local record is
mapped with them; the fit residual and the helicity off-diagonal of the mapped records are printed (the latter is a check of the C3H
fix: a wrong o_H sign would not be helicity-diagonal).
Variants (argv[1]):
  n1val : + local n1 at x = 0.85, 0.87, 0.90, 0.92, 0.93                   (Gamma reference for the n3 comparison)
  n3val : + local patched n3 at the same nodes                            (compare with n1val, node by node, in c2)
  ovl   : batch-3 n1 at x = 0.9198, 0.9276 REPLACED by local patched n3   (n3 vs external n1 at the same nu)
  edge  : + the patched-n3 GK edge nodes (x = 0.9357 ... 0.9746)          (the k = 2 AND k = 3 edge filled; single seed, M = 10)
  p1/p3 : + x = 0.92 at higher precision: n1 M = 13 (p1) / patched n3 M = 12, independent seed 7000 (p3)
Usage: python3 nf_gamma.py VARIANT  -> stdout, cache/gamma/<variant>/gamma_J26.{json,txt}
"""
import nf_guard  # noqa: F401  (first)
import os, sys, pickle, glob
sys.modules['e2_guard'] = sys.modules['er_guard'] = sys.modules['nf_guard']        # e2's own guard (root edger2/cache) not installed
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for d in ('edger2', 'edger', 'highspin'):
    sys.path.insert(0, os.path.join(CODE, d))
import warnings; warnings.filterwarnings('ignore')
import numpy as np
VAR = sys.argv[1]
J, D, DREF = 26, 1026.0, 20000.0
OUTD = os.path.join(HERE, 'cache', 'gamma', VAR)
os.makedirs(OUTD, exist_ok=True)
import e2_gamma as G
G.CACHE = OUTD
G.FIX = os.path.join(OUTD, 'ptsfix')
EB = G.EB
ANCH = (942.006243, 877.94542)
SRC = {'n1val': ['val_n1'], 'n3val': ['val_n3'], 'ovl': ['ovl'], 'edge': ['edge', 'edge2', 'edge3'], 'p1': ['prec_n1'], 'p3': ['prec_n3']}[VAR]


def lp(sub, kind):
    f = os.path.join(HERE, 'cache', sub, 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D))
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def full(r):
    return np.exp(r['logN'])*np.asarray(r['reps'])            # (R, 2, 2) full-integral replicas


def fitAB(X, Y, it=200):
    """min sum_i |A X_i B - Y_i|^2 (alternating LS); X, Y: (n, 2, 2)."""
    s = np.abs(Y).max()/np.abs(X).max(); X = X*s
    B = np.eye(2, dtype=complex)
    for _ in range(it):
        Xb = X @ B
        A = np.einsum('nij,nkj->ik', Y, Xb.conj()) @ np.linalg.inv(np.einsum('nij,nkj->ik', Xb, Xb.conj()))
        Z = A @ X
        B = np.linalg.inv(np.einsum('nji,njk->ik', Z.conj(), Z)) @ np.einsum('nji,njk->ik', Z.conj(), Y)
    res = np.linalg.norm(A @ X @ B - Y)/np.linalg.norm(Y)
    return A*s, B, res


LOG = []
_fix0 = EB.fix


def fix(J_, D_, HN, Tloc, perm=(None, None), src=None, dst=None):
    rep = _fix0(J_, D_, HN, Tloc, perm, src, dst)
    P = {k: pickle.load(open(os.path.join(dst, 'hc_pts_%s_J%d_D%g.pkl' % (k, J_, D_)), 'rb')) for k in ('EE', 'C3H')}
    A_ = {k: lp('anchor', k) for k in ('EE', 'C3H')}
    LOG.clear()
    for DD in (D_, DREF):
        X, Y = [], []
        for kind in ('EE', 'C3H'):
            for nu in ANCH:
                key = (DD, round(nu, 9))
                fl, fe = full(A_[kind][key]), full(P[kind][key])
                X += list(fl) + [fl.mean(0)]; Y += list(fe) + [fe.mean(0)]
        Am, Bm, res = fitAB(np.array(X), np.array(Y))
        (b, logK, _), _ = HN.shadow(DD, J_, HN.RESID.get((float(DD), J_, 'EE')))
        T = Tloc[DD]; Ti = np.linalg.inv(T); Bi = np.linalg.inv((np.asarray(b) @ T).T)
        LOG.append(f"  Delta={DD:g}: anchor fit residual {res:.1e}")
        for sub in SRC:
            for kind in ('EE', 'C3H'):
                for key, r in lp(sub, kind).items():
                    if key[0] != DD:
                        continue
                    f2 = np.array([Am @ q @ Bm for q in full(r)])
                    q = dict(r); q['reps'] = f2*np.exp(-r['logN']); q['E'] = q['reps'].mean(0)
                    q['err'] = q['reps'].std(axis=0, ddof=1)/np.sqrt(len(q['reps']))
                    C = Ti @ q['E'] @ Bi
                    off = np.hypot(abs(C[0, 1]), abs(C[1, 0]))/np.linalg.norm(C)
                    LOG.append(f"    {sub} {kind} Delta={DD:g} x={key[1]/G.HN.mphys(D_, J_):.4f}: injected{' (replaces batch-3 n1)' if key in P[kind] else ''}; "
                               f"helicity off-diagonal / ||C|| = {off:.1e}")
                    P[kind][key] = q
    for k in ('EE', 'C3H'):
        pickle.dump(P[k], open(os.path.join(dst, 'hc_pts_%s_J%d_D%g.pkl' % (k, J_, D_)), 'wb'))
    return rep


EB.fix = fix
res = G.main(J)
print('\n'.join(LOG))
