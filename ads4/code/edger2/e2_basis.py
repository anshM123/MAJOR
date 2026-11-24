"""
edger2/e2_basis.py (EDGE_R2.md §1.1): STRUCTURE-BASIS REPAIR of the external even-J records.

Finding: the external batch-3 even-J records are NOT in the conserved-structure basis of this machine.  l5_spin.conserved returns an
SVD basis of the <TTO> null space (hc_points' docstring: the points "must be in the SAME structure basis as the existing K2 rows they
are injected into"); on the external machine that basis differs from the local one by a general 2x2 transformation (not only column
signs: e2_report §1), so hc_nodes.to_record with the local shadow b gives matrices that are NOT diagonal in the chain's helicity basis
(off/diag ~ 0.1-2 instead of <~ 1e-5 as for the local J = 16 row), and the chain's Gamma_{k,+-} is garbage (x10-1000).
(Odd J is 1x1: a basis change is a constant factor and cancels in a_h -- the odd external rows of EXT_FOLD and J = 33 here are fine.)

Repair (exact in principle, no new physics input): in ANY structure basis the all-orders matrix of a state is E(nu) = U C(nu) V^T with
C diagonal in helicity and U, V constant (the helicity decomposition is nu-independent; locally the direct records are diagonal in the
K2 basis to <~ 5e-5).  U, V are found from the external records themselves by joint diagonalisation over all nodes (EE kind; angles of
the columns of U and V, least squares on |offdiag|/sqrt|diag diag|), separately for the state row and for Delta_ref = 2e4 (whose nodes
x_ref <= 0.05 straddle the reference zero of D_++, which separates the two helicities).  The records are then mapped into the LOCAL
basis, E' = T U^-1 E V^-T (b T)^T (T = the chain's l8 helicity basis, b = the local shadow), so that the chain runs unchanged:
Ti M0' Ti^T = U^-1 E V^-T.  The per-helicity constants of U, V are arbitrary; they are absorbed in the chain's a_h (common to k = 2, 3:
the same U, V are applied to the C3H records).  The helicity LABEL of each column is fixed afterwards (e2_gamma: the permutation for
which both Gamma_h are smallest on x < 0.6; the alternative is reported).
Checks reported: residual off/diag over all nodes, for EE (fitted) and C3H (not fitted).
Usage: imported by e2_gamma (fix(J, D, perm_state, perm_ref) -> directory of repaired point files)
"""
import e2_guard  # noqa: F401
import os, pickle
import numpy as np
from scipy.optimize import least_squares
HERE = os.path.dirname(os.path.abspath(__file__))
DREF = 20000.0


def _UV(p):
    a, b, c, d = p
    U = np.array([[np.cos(a), np.cos(b)], [np.sin(a), np.sin(b)]])
    V = np.array([[np.cos(c), np.cos(d)], [np.sin(c), np.sin(d)]])
    return U, V


def _diag(E, U, V):
    return np.linalg.solve(U, E) @ np.linalg.inv(V).T


def offrel(Es, U, V):
    """max |off-diagonal| / ||C|| per node (C = U^-1 E V^-T)."""
    out = []
    for E in Es:
        C = _diag(E, U, V)
        out.append(max(abs(C[0, 1]), abs(C[1, 0]))/(np.linalg.norm(C) + 1e-300))
    return np.array(out)


def joint(Es):
    """U, V (unit columns) with U^-1 E_i V^-T ~ diagonal for all i."""
    Es = [E/np.abs(E).max() for E in Es]
    best = None
    n = len(Es)
    for i in range(0, n, max(1, n//6)):
        for j in range(n - 1, i, -max(1, n//6)):
            try:
                lam, Vt = np.linalg.eig(np.linalg.solve(Es[i], Es[j]))
                lam2, U = np.linalg.eig(Es[j] @ np.linalg.inv(Es[i]))
            except np.linalg.LinAlgError:
                continue
            if abs(lam[0] - lam[1]) < 1e-6*abs(lam).max():
                continue
            o = [int(np.argmin(abs(lam2 - l))) for l in lam]
            U = U[:, o]
            V = np.linalg.inv(Vt).T
            ang = []
            for M in (U, V):
                for k in range(2):
                    v = M[:, k]; v = v/v[np.argmax(abs(v))]; v = v.real
                    ang.append(np.arctan2(v[1], v[0]))
            p0 = np.array(ang)
            r = offrel(Es, *_UV(p0))
            sc = float(np.max(r))
            if best is None or sc < best[0]:
                best = (sc, p0)
    p0 = best[1]

    def res(p):
        U, V = _UV(p); out = []
        for E in Es:
            C = _diag(E, U, V); s = np.linalg.norm(C) + 1e-300
            out += [abs(C[0, 1])/s, abs(C[1, 0])/s]
        return np.array(out)
    sol = least_squares(res, p0, xtol=1e-15, ftol=1e-15, gtol=1e-15, loss='soft_l1', f_scale=1e-5)
    U, V = _UV(sol.x)
    return U, V


def fix(J, D, HN, Tloc, perm=(None, None), src=None, dst=None):
    """repaired copies of the (EE, C3H) point files of state D: records of D and of DREF mapped into the local basis.
    Tloc: {Delta: T (columns ++, +-)} (the chain's l8 basis); perm: per Delta (state, ref) a column permutation (0,1) or (1,0)."""
    src = src or HN.C
    os.makedirs(dst, exist_ok=True)
    P = {k: pickle.load(open(os.path.join(src, 'hc_pts_%s_J%d_D%g.pkl' % (k, J, D)), 'rb')) for k in ('EE', 'C3H')}
    rep = {}
    for i, DD in enumerate((D, DREF)):
        keys = sorted(k for k in P['EE'] if k[0] == DD and k in P['C3H'])
        Es = [P['EE'][k]['E'] for k in keys]
        U, V = joint(Es)
        pm = perm[i] or (0, 1)
        U, V = U[:, list(pm)], V[:, list(pm)]
        (b, logK, res_), _ = HN.shadow(DD, J, HN.RESID.get((float(DD), J, 'EE')))
        T = Tloc[DD]
        A = T @ np.linalg.inv(U)
        Bm = np.linalg.inv(V).T @ (np.asarray(b) @ T).T
        rep[DD] = dict(off_EE=offrel(Es, U, V), off_C3H=offrel([P['C3H'][k]['E'] for k in keys], U, V),
                       x=np.array([k[1] for k in keys]), U=U, V=V)
        for kind in ('EE', 'C3H'):
            for k in keys:
                r = dict(P[kind][k])
                r['E'] = A @ r['E'] @ Bm
                r['reps'] = np.array([A @ q @ Bm for q in r['reps']])
                r['err'] = np.abs(A) @ np.asarray(r['err']) @ np.abs(Bm)     # (only for the failure screen: conservative)
                P[kind][k] = r
    for kind in ('EE', 'C3H'):
        pickle.dump(P[kind], open(os.path.join(dst, 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)), 'wb'))
    return rep
