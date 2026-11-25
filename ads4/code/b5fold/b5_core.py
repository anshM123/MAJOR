"""
b5_core (BATCH5_FOLD.md): shared set-up for folding the batch-5 external direct points into the Delta_g = 1000 near-gap budgets.
A copy of extfold/ef_core.py (merged point reader for smallnu/sn_core -> highspin/hs_core -> neargap2/ng_val -> spinfunc chain) with
  * the point sets moved to b5fold/cache/pts/<set>/ (b5_merge.py),
  * every output redirected to b5fold/cache (sn_core.SNC / HSC, hs_core.HSC),
  * the K2 EE state records of the even own rows (b5_k2.py; for the chain's l8 helicity basis) added to the row loader,
  * the even <TTO> shadows of the own rows (1000 + J, J), computed here (b5_shadow.py, cache/b5_shadow_k2.pkl), appended to
    hc_nodes' even candidate list (no K2 row exists at Delta = 1000 + J for even J = 18 ... 32).
Every project module is imported read-only; b5_guard forbids any write outside b5fold/cache.
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import b5_guard  # noqa: F401  (FIRST: stands in for ng_guard / hs_guard / sn_guard / ef_guard / nf_guard)
SN = os.path.join(HERE, '..', 'smallnu')
HSD = os.path.join(HERE, '..', 'highspin')
sys.path.insert(0, HSD); sys.path.insert(0, SN); sys.path.insert(0, HERE)
import numpy as np
import sn_core as S                    # read-only import (it patches hc_nodes.pts; patched again below)
H = S.H
HN = H.HN
V = H.V
B5C = os.path.join(HERE, 'cache')
PTS = os.path.join(B5C, 'pts')
SET = ['all']
HS_CACHE = H.HSC                       # highspin/cache (read-only)
NUWIN = 80.0
SHF = os.path.join(B5C, 'b5_shadow_k2.pkl')


def pfile(kind, J, D, st=None):
    return os.path.join(PTS, st or SET[0], 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D))


_CACHE = {}


def _pts(kind, J, D):
    st = 'base' if S.BASE[0] else SET[0]
    f = pfile(kind, J, D, st)
    if not os.path.exists(f):
        return {}
    k = (f, os.path.getmtime(f))
    if k not in _CACHE:
        _CACHE[k] = pickle.load(open(f, 'rb'))
    return dict(_CACHE[k])


def sn_nus(J, D, sub=None):
    P = _pts('EE', J, D)
    return sorted(k[1] for k in P if k[0] == D and (H.DREF, k[1]) in P and k[1] < NUWIN)


def rel_err(J, D, nu, sub=None):
    e = 0.0
    for kind in ('EE', 'C3H'):
        P = _pts(kind, J, D)
        for DD in (D, H.DREF):
            for k, r in P.items():
                if k[0] == DD and abs(k[1] - nu) < 1e-4*max(nu, 1.0):
                    e = max(e, float(np.max(np.abs(r['err'])/np.maximum(np.abs(r['E']), 1e-300))))
    return e


HN.pts = _pts
S._pts = _pts
S.sn_nus = sn_nus
S.rel_err = rel_err
S.SNC = B5C
S.HSC = B5C
H.HSC = B5C

_cands0 = HN._cands


def _cands(odd):
    d = _cands0(odd)
    if odd:
        return d
    if 'even_b5' not in HN._SC:
        dd = {k: list(v) for k, v in d.items()}
        if os.path.exists(SHF):
            for k, v in pickle.load(open(SHF, 'rb')).items():
                dd.setdefault(k, []).append(('b5_shadow_k2.pkl', v))
        HN._SC['even_b5'] = dd
    return HN._SC['even_b5']


HN._cands = _cands

# the K2 EE state records of the even own rows (b5_k2.py: l8 helicity basis only, not row nodes) are read by the chain's row loader
import glob as _glob
_W = H.NC.W
_pk0 = _W._pk


def _pk(pats):
    R = _pk0(pats)
    for p in pats:
        for fn in sorted(_glob.glob(os.path.join(B5C, p))):
            R.update(pickle.load(open(fn, 'rb')))
    return R


_W._pk = _pk
H.NC._pk = _pk


def mphys(D, J):
    return HN.mphys(D, J)


def twist(J, Dg=1000.0):
    return mphys(Dg + J, J)/H.E.lam_of(Dg)


# The chain's Hampel despike (w2_chain._grid, W2_DESPIKE) is meant for K2 rows: injected DIRECT rows are exempt (the HC protocol,
# HC_DESPIKE_INJ=0; neargap2/ng_val._clean's docstring).  ng_bind.run_state builds its grids through Q.even_grid -> Q.A.grid directly,
# so there the filter also hit the direct own row (NEARGAP_FILL's J = 14: 'despiked 1014 [0.925]'; here up to 5 edge nodes per row,
# e.g. J = 26: 0.8226, 0.8362, 0.923, 0.925, 0.9387).  Default here: the rows being injected (V.INJ[J]) are exempt, every K2 row is
# despiked exactly as before (w2_chain's rule, copied).  B5_HAMPEL_DIRECT=1 restores the as-run behaviour (comparison).
HAMPEL_DIRECT = os.environ.get('B5_HAMPEL_DIRECT', '0') == '1'
_Q = H.NC.Q
_Ag = _Q.A.grid


def _hampel(g, J, skip):
    for D, r in g.items():
        if D in skip:
            continue
        n = len(r['xs']); bad = np.zeros(n, bool)
        for k in [kk for kk in r if isinstance(kk, tuple) and len(kk) == 2 and kk[1] in ('++', '+-')]:
            v = np.asarray(r[k], float)
            scale = np.median(np.abs(v)) + 1e-12
            for i in range(n):
                w = np.r_[v[max(0, i - 3):i], v[i + 1:i + 4]]
                med = np.median(w); mad = np.median(np.abs(w - med))
                if abs(v[i] - med) > 6*max(mad, scale):
                    bad[i] = True
        bad[0] = False
        if bad.any():
            _W.DESPIKE_LOG.append((J, D, [round(float(x), 4) for x in r['xs'][bad]]))
            for k, v in list(r.items()):
                if isinstance(v, np.ndarray) and v.shape == bad.shape:
                    r[k] = v[~bad]
    return g


def _grid_exempt(J, K, *a, **k):
    inj = {float(D) for D in V.INJ.get(J, [])}
    if HAMPEL_DIRECT or not inj or os.environ.get('W2_DESPIKE', '1') != '1':
        return _Ag(J, K, *a, **k)
    os.environ['W2_DESPIKE'] = '0'
    try:
        g = _Ag(J, K, *a, **k)
    finally:
        os.environ['W2_DESPIKE'] = '1'
    return _hampel(g, J, inj)


_Q.A.grid = _grid_exempt
