"""
fw_fix (FOLD_WINDOW_LARGEJ.md §1.6): load-time consistency of the structure basis across the fw records.
The conserved <TTO> null space (l5_spin.conserved) is returned in an SVD-arbitrary orthonormal basis that depends on the random sample
points; fw_frame.patch_conserved(_clean) (z normalisation, needed at J >= 150) changes the points and hence rotates that basis.  Records
computed before the patch was wired in (10:10 UTC; 'unpatched') and after ('patched') therefore differ by fixed 2x2 maps, and K2 records
store M = S b^-T with whatever shadow b was current.  Every row must be in ONE basis (the generalised helicity basis absorbs one fixed map
per side and row, not a mixture).  Fixes (exact linear algebra, nothing on disk changed):
  * K2 records: the b actually used is identified by its stored shadow_resid among the known shadow caches; M_o <- M_o b_used^T b_cur^-T
    with b_cur = cache/fw_shadowk2.pkl (all computed in the patched basis; the first-run shadows are kept in cache/fw_shadowrun1.pkl;
    neither file name matches other agents' '*shadow_k2*' globs);
  * direct points listed in UNPATCHED (computed before the patch): E <- P^-T E Pt^-1, P = pinv(sol_pat) sol_unp (Delta side),
    Pt likewise (3 - Delta side), for E, err (abs of transformed err, approximate), reps.
"""
import os, sys, glob, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
C = os.path.join(HERE, 'cache')
UNPATCHED = {(100, 400.0, 16.7)}          # (J, Delta, nu) of direct points with the unpatched structure basis.  The (400, 100) run's
# reference at (2e4, 16.7) was ALSO unpatched, but it was overwritten at ~18:30 UTC by the (1100, 100) run's (patched) reference at the
# same key (FOLD_WINDOW_LARGEJ §2.6), so it is no longer converted.
_B = {}


def _cands():
    if not _B:
        fs = sorted(glob.glob(os.path.join(C, '*shadow_k2*.pkl'))) + [os.path.join(C, f) for f in ('fw_shadowrun1.pkl', 'fw_shadowk2.pkl')
                                                                     if os.path.exists(os.path.join(C, f))]
        for f in fs:
            for k, v in pickle.load(open(f, 'rb')).items():
                _B.setdefault(k, []).append((os.path.basename(f), v))
    return _B


def bcur(D, J):
    return pickle.load(open(os.path.join(C, 'fw_shadowk2.pkl'), 'rb'))[(float(D), int(J))]


def fix_k2(R, J, log=None):
    out = {}
    cur = {}
    for (D, nu), r in R.items():
        if D not in cur:
            cur[D] = bcur(D, J)
        bc = cur[D]
        res = np.asarray(r.get('shadow_resid', []), float).ravel()
        if res.size and np.allclose(res, np.asarray(bc[2], float).ravel(), rtol=1e-12, atol=0):
            out[(D, nu)] = r; continue
        used = [v for f, v in _cands().get((float(D), int(J)), []) if np.allclose(res, np.asarray(v[2], float).ravel(), rtol=1e-12, atol=0)]
        if not used:
            raise RuntimeError(f"fw_fix: no shadow matches record ({D}, {nu}) J={J}")
        bu = used[0][0]
        X = bu.T @ np.linalg.inv(bc[0]).T
        r2 = dict(r)
        for o in ('M0', 'M2', 'M4', 'M6'):
            if o in r2:
                r2[o] = r2[o] @ X
        r2['logK'] = r2['logK'] + (bc[1] - used[0][1])
        r2['_fwfix'] = 'b'
        out[(D, nu)] = r2
        if log is not None:
            log.append((J, D, round(nu, 3)))
    return out


_SOL = {}


def _sols(D, J):
    """(sol_unpatched, sol_patched) for dimension D (l5_spin.conserved post-processing applied to both)."""
    key = (float(D), int(J))
    if key not in _SOL:
        import l5_spin as S
        import fw_frame as FW
        L = S._l4r()
        def post(r):
            sol = r['sol']
            U, sv_, _ = np.linalg.svd(np.hstack([sol.real, sol.imag]), full_matrices=False)
            sol = U[:, :int(np.sum(sv_ > 1e-8*sv_[0]))]
            return sol/np.abs(sol).max(axis=0, keepdims=True)
        orig_pt, orig_ce = getattr(L, '_fw_orig_rp', None), getattr(L, '_fw_orig_ce', None)
        if getattr(L, '_fw_norm', False):
            raise RuntimeError('fw_fix._sols must run before fw_frame.patch_conserved in this process')
        su = post(L.conserved_even(float(D), int(J)))
        FW.patch_conserved_clean(S)
        sp = post(L.conserved_even(float(D), int(J)))
        _SOL[key] = (su, sp)
    return _SOL[key]


def precompute_sols(pairs):
    """call once, before any fw_frame patch, for all (D, J) needed (both D and 3 - D)."""
    import l5_spin as S
    L = S._l4r()
    import fw_frame as FW
    out = {}
    for D, J in pairs:
        for DD in (D, 3.0 - D):
            def post(r):
                sol = r['sol']
                U, sv_, _ = np.linalg.svd(np.hstack([sol.real, sol.imag]), full_matrices=False)
                sol = U[:, :int(np.sum(sv_ > 1e-8*sv_[0]))]
                return sol/np.abs(sol).max(axis=0, keepdims=True)
            out[(float(DD), int(J))] = [post(L.conserved_even(float(DD), int(J)))]
    FW.patch_conserved_clean(S)
    for k in out:
        out[k].append(post(L.conserved_even(*k)))
        _SOL[k] = tuple(out[k])
    pickle.dump(_SOL, open(os.path.join(C, 'fw_sols.pkl'), 'wb'))
    return _SOL


def load_sols():
    if not _SOL:
        f = os.path.join(C, 'fw_sols.pkl')
        if os.path.exists(f):
            _SOL.update(pickle.load(open(f, 'rb')))
    return _SOL


def fix_direct(r, J, D, nu):
    if (int(J), float(D), round(float(nu), 1)) not in UNPATCHED:
        return r
    S_ = load_sols()
    su, sp = S_[(float(D), int(J))]; tu, tp = S_[(3.0 - float(D), int(J))]
    P = np.linalg.pinv(sp) @ su; Pt = np.linalg.pinv(tp) @ tu
    resid = max(np.abs(sp @ P - su).max(), np.abs(tp @ Pt - tu).max())
    assert resid < 1e-8, resid
    A, B = np.linalg.inv(P).T, np.linalg.inv(Pt)
    r2 = dict(r)
    r2['E'] = A @ r['E'] @ B
    r2['reps'] = np.array([A @ x @ B for x in r['reps']])
    r2['err'] = np.abs(A) @ r['err'] @ np.abs(B)
    r2['_fwfix'] = 'sol'
    return r2


if __name__ == '__main__':
    pairs = [(400.0, 100), (20000.0, 100)]
    S_ = precompute_sols(pairs)
    for k, (u, p) in S_.items():
        P = np.linalg.pinv(p) @ u
        print(k, 'P =', np.round(P, 6).ravel(), 'resid', np.abs(p @ P - u).max())
