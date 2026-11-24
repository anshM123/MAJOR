"""
b5_core (BATCH5_FOLD.md): shared set-up for folding the batch-5 external direct points into the Delta_g = 1000 near-gap budgets.
A copy of extfold/ef_core.py (merged point reader for smallnu/sn_core -> highspin/hs_core -> neargap2/ng_val -> spinfunc chain) with
  * the point sets moved to b5fold/cache/pts/<set>/ (b5_merge.py),
  * every output redirected to b5fold/cache (sn_core.SNC / HSC, hs_core.HSC),
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


def mphys(D, J):
    return HN.mphys(D, J)


def twist(J, Dg=1000.0):
    return mphys(Dg + J, J)/H.E.lam_of(Dg)
