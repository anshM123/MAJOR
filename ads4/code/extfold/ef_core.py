"""
ef_core (EXT_FOLD.md): fold the external direct points (ads4/ext, jobs A and B) into the near-gap positivity budgets at Delta_g = 1000,
through the existing analysis code imported READ-ONLY: smallnu/sn_core (-> highspin/hs_core -> neargap2/ng_val -> spinfunc chain),
highspin/hs_odd, smallnu/sn_odd, smallnu/sn_window.

Point sources (all read-only; merged by ef_merge.py into extfold/cache/pts/<set>/hc_pts_<kind>_J<J>_D<D>.pkl):
  set 'all'  : every direct point of the state (spinfunc/cache [HC], neargap2/cache, highspin/cache, smallnu/cache top level, and the
               ext job pickles; for J = 29 the n3 small-nu records POOLED over all distinct RQMC streams, ef_merge.pool_n3)
  set 'base' : the points of the HIGHSPIN baseline grid (highspin/cache only) -- the normalisation baseline of sn_core (BASE[0] = True)
  set 'old'  : 'all' without the ext points (J = 29: smallnu's seed-4242 n3 records unpooled) -- for the old -> new comparison
SET[0] selects the set used when BASE[0] is False.  Outputs go to extfold/cache (ef_guard forbids any other write).
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ef_guard  # noqa: F401  (FIRST: stands in for ng_guard / hs_guard / sn_guard)
SN = os.path.join(HERE, '..', 'smallnu')
HSD = os.path.join(HERE, '..', 'highspin')
sys.path.insert(0, HSD); sys.path.insert(0, SN); sys.path.insert(0, HERE)
import numpy as np
import sn_core as S                    # read-only import (it patches hc_nodes.pts; we patch it again below)
H = S.H
HN = H.HN
EFC = os.path.join(HERE, 'cache')
PTS = os.path.join(EFC, 'pts')
SET = ['all']
HS_CACHE = H.HSC                       # highspin/cache (read-only; HIGHSPIN's result files)
NUWIN = 80.0                           # direct nodes below this nu are 'window' (sn) nodes for sn_core's labels / rel-err filter


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
S.SNC = EFC                            # sn_odd writes S.SNC/sn_odd_J*.{txt,json} and reads S.SNC/hs_odd_J*.json (FULL normalisation)
S.HSC = EFC
H.HSC = EFC                            # hs_odd (run via runpy) writes H.HSC/hs_odd_J*.{txt,json}

# odd shadows computed here (ef_shadow.py: J = 7, Delta = 1007, which has no K2 row) appended to hc_nodes' candidate list
_cands0 = HN._cands


def _cands(odd):
    d = _cands0(odd)
    if not odd:
        return d
    if 'odd_ef' not in HN._SC:
        dd = {k: list(v) for k, v in d.items()}
        f = os.path.join(EFC, 'ef_p2shadow.pkl')
        if os.path.exists(f):
            for k, v in pickle.load(open(f, 'rb')).items():
                dd.setdefault(k, []).append(('ef_p2shadow.pkl', v))
        HN._SC['odd_ef'] = dd
    return HN._SC['odd_ef']


HN._cands = _cands
