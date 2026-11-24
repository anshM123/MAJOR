"""
nf_core (NEARGAP_FILL.md): shared set-up for the Delta_g = 1000 near-gap fill (J = 9, 10, 11, 13, 14, 15).
Everything is the unchanged project machinery, imported read-only:
  highspin/hs_core -> neargap2/ng_val -> neargap2/ng_chain -> spinfunc (w2/q3 chain in EE mode, ee1000 Jacobi fix).
Redirections (process-local; nf_guard forbids every write outside ngfill/cache):
  * row loader (w2_chain._pk): spinfunc/cache + neargap2/cache + ngfill/cache (the new own-mass K2 rows of nf_rows.py);
  * direct points (hc_nodes.C): ngfill/cache (nf_points.py);
  * shadow candidates (hc_nodes._cands via ng_val): also ngfill/cache's nf_p2shadow.pkl / nf_shadow_k2.pkl.
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import nf_guard  # noqa: F401  (first: registers itself as ng_guard / hs_guard / ...)
HS = os.path.join(HERE, '..', 'highspin')
NG = os.path.join(HERE, '..', 'neargap2')
SF = os.path.join(HERE, '..', 'spinfunc')
NFC = os.path.join(HERE, 'cache')
sys.path.insert(0, SF); sys.path.insert(0, NG); sys.path.insert(0, HS)
import numpy as np
import hs_core as H                      # read-only
V, NC, HN = H.V, H.NC, H.HN
E, Q, W = NC.E, NC.Q, NC.W
NGC = os.path.join(NG, 'cache')
SFC = NC.SFC


def _pk(pats):
    R = {}
    for p in pats:
        fns = sorted(glob.glob(os.path.join(SFC, p)))
        if p.startswith('w1_p3_'):
            fns += sorted(glob.glob(os.path.join(SFC, 'ee1000_p3_' + p[len('w1_p3_'):])))
        if not NC.NOOWN:
            fns += sorted(glob.glob(os.path.join(NGC, p))) + sorted(glob.glob(os.path.join(NFC, p)))
        for fn in fns:
            R.update(pickle.load(open(fn, 'rb')))
    return R


W._pk = _pk
NC._pk = _pk
HN.C = NFC
HN._SC.clear()
_cands_prev = HN._cands                  # ng_val's wrapper (project + neargap2 shadows)


def _cands(odd):
    d = _cands_prev(odd)
    key = 'odd_nf' if odd else 'even_nf'
    if key not in HN._SC:
        dd = {k: list(v) for k, v in d.items()}
        f = os.path.join(NFC, 'nf_p2shadow.pkl' if odd else 'nf_shadow_k2.pkl')
        if os.path.exists(f):
            for k, v in pickle.load(open(f, 'rb')).items():
                dd.setdefault(k, []).append((os.path.basename(f), v))
        HN._SC[key] = dd
    return HN._SC[key]


HN._cands = _cands


def mphys(D, J):
    return HN.mphys(D, J)


def twist(J, Dg=1000.0):
    return mphys(Dg + J, J)/E.lam_of(Dg)
