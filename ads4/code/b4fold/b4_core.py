"""
b4_core (BATCH4_FOLD.md): shared set-up for folding the batch-4 external direct points into the Delta_g = 2000 near-gap budgets.
A copy of b5fold/b5_core.py (merged point reader for smallnu/sn_core -> highspin/hs_core -> neargap2/ng_val -> spinfunc chain) with
  * the point sets moved to b4fold/cache/pts/<set>/ (b4_merge.py; 'fix' = even J mapped to the local structure basis, b4_basis.py),
  * every output redirected to b4fold/cache (sn_core.SNC / HSC, hs_core.HSC),
  * Delta_g = 2000 (GAP_NEARGAP's in-process redirections of gapnear/gn_bind.py, re-typed):
      - the row loader w2_chain._pk also reads gapnear/cache/rows (the own-mass K2 rows 2000 + J, Delta_ref = 2e4: the reference the
        batch-4 direct points use) unless the chain is asked for the published rows only (NG_NOOWN), and the K2 EE state records of the
        seam row 2020 (b4_k2.py; for the even chain's l8 helicity basis, not row nodes);
      - ng_bind.gauge_factor at Delta_g = 2000 uses the min-norm weight gapscale/cache/weights_2000.npz (read only), as gn_bind;
  * the shadows of the 2000 rows: gapnear's gn_shadow_k2 / gn_p2shadow and this note's (2020, 20), (2023, 23, 6) (b4_shadow.py)
    appended to hc_nodes' candidate lists,
  * the Hampel exemption for injected DIRECT rows (BATCH5_FOLD §2.2 item 2), copied unchanged.
Every project module is imported read-only; b4_guard forbids any write outside b4fold/cache.
"""
import os, sys, pickle, glob as _glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import b4_guard  # noqa: F401  (FIRST: stands in for ng_guard / hs_guard / sn_guard / ef_guard / nf_guard / b5_guard)
SN = os.path.join(HERE, '..', 'smallnu')
HSD = os.path.join(HERE, '..', 'highspin')
sys.path.insert(0, HSD); sys.path.insert(0, SN); sys.path.insert(0, HERE)
import numpy as np
import sn_core as S                    # read-only import (it patches hc_nodes.pts; patched again below)
H = S.H
HN = H.HN
V = H.V
B4C = os.path.join(HERE, 'cache')
PTS = os.path.join(B4C, 'pts')
SET = ['fix']
GNR = os.path.abspath(os.path.join(HERE, '..', 'gapnear', 'cache', 'rows'))
GSW = os.path.abspath(os.path.join(HERE, '..', 'gapscale', 'cache'))
NUWIN = 80.0
DG = 2000.0
S.FULL[0] = True                       # no baseline grid: every row's normalisation is that of the grid with all its direct points


def pfile(kind, J, D, st=None):
    return os.path.join(PTS, st or SET[0], 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D))


_CACHE = {}


def _pts(kind, J, D):
    f = pfile(kind, J, D)
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
S.SNC = B4C
S.HSC = B4C
H.HSC = B4C

# ------------------------------------------------------------------------------------------------ shadows of the 2000 rows
_cands0 = HN._cands


def _cands(odd):
    d = _cands0(odd)
    key = 'odd_b4' if odd else 'even_b4'
    if key not in HN._SC:
        dd = {k: list(v) for k, v in d.items()}
        for f in ((os.path.join(GNR, 'gn_p2shadow.pkl'), os.path.join(B4C, 'b4_p2shadow.pkl')) if odd else
                  (os.path.join(GNR, 'gn_shadow_k2.pkl'), os.path.join(B4C, 'b4_shadow_k2.pkl'))):
            if os.path.exists(f):
                for k, v in pickle.load(open(f, 'rb')).items():
                    dd.setdefault(k, []).append((os.path.basename(f), v))
        HN._SC[key] = dd
    return HN._SC[key]


HN._cands = _cands

# ------------------------------------------------------------------------------------------------ row loader: gapnear 2000 rows + K2 EE state records
NC = H.NC
_W = NC.W
_pk0 = _W._pk


def _pk(pats):
    R = _pk0(pats)
    for p in pats:
        fns = sorted(_glob.glob(os.path.join(B4C, p)))
        if not NC.NOOWN:
            fns = sorted(_glob.glob(os.path.join(GNR, p))) + fns
        for fn in fns:
            R.update(pickle.load(open(fn, 'rb')))
    return R


_W._pk = _pk
NC._pk = _pk

# ------------------------------------------------------------------------------------------------ gauge factor at Delta_g = 2000
import ng_bind as NB                   # read-only import
_WTS = {2000: dict(np.load(os.path.join(GSW, 'weights_2000.npz')))}
_gf_orig = NB.gauge_factor


def gauge_factor(Dg, nu):
    key = int(round(Dg))
    if key not in _WTS:
        return _gf_orig(Dg, nu)
    w = _WTS[key]; g = w['nu']; X = w['X'].real; phi = w['phi']
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.abs(X/phi)
    r = np.where(np.isfinite(r), r, 1.0)
    out = np.ones_like(nu)
    for i, v in enumerate(nu):
        if v <= 25.0:
            sel = np.abs(g - v) <= 0.5
            out[i] = max(1.0, float(r[sel].max()) if sel.any() else 1.0)
    return out


NB.gauge_factor = gauge_factor


def mphys(D, J):
    return HN.mphys(D, J)


def twist(J, Dg=DG):
    return mphys(Dg + J, J)/H.E.lam_of(Dg)


# ------------------------------------------------------------------------------------------------ Hampel exemption (b5_core, copied)
HAMPEL_DIRECT = os.environ.get('B4_HAMPEL_DIRECT', '0') == '1'
_Q = NC.Q
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


DROPX = {}          # (J, D) -> x of the own row's nodes that are not used (failed / reference-zero rule; set by b4_state)
DROP_LOG = []


def _drop_own(g, J, inj):
    """On an injected own row that also has a K2 row (J = 12, 16 at 2000 + J): a direct node removed by the RQMC-failure or reference-
    zero rule must not fall back to the K2 node at the same x (that K2 node is the one the rule is about; the K2-only chain removed it
    too, by its Hampel filter or refbad), and a node with a non-finite entry is dropped (ng_val._clean's rule)."""
    for D in inj:
        if D not in g:
            continue
        r = g[D]; n = len(r['xs'])
        xb = np.array(sorted(DROPX.get((J, D), ())), float)
        bad = np.array([len(xb) > 0 and np.min(np.abs(xb - x)) < 1e-5 for x in r['xs']])
        for k, v in r.items():
            if isinstance(v, np.ndarray) and v.shape[:1] == (n,) and np.issubdtype(v.dtype, np.number):
                bad |= ~np.isfinite(np.asarray(v, complex).reshape(n, -1)).all(1)
        if bad.any():
            DROP_LOG.append((J, D, [round(float(x), 4) for x in r['xs'][bad]]))
            g[D] = {kk: (v[~bad] if isinstance(v, np.ndarray) and v.shape[:1] == (n,) else v) for kk, v in r.items()}
    return g


def _grid_exempt(J, K, *a, **k):
    inj = {float(D) for D in V.INJ.get(J, [])}
    if HAMPEL_DIRECT or not inj or os.environ.get('W2_DESPIKE', '1') != '1':
        return _drop_own(_Ag(J, K, *a, **k), J, inj)
    os.environ['W2_DESPIKE'] = '0'
    try:
        g = _Ag(J, K, *a, **k)
    finally:
        os.environ['W2_DESPIKE'] = '1'
    return _drop_own(_hampel(g, J, inj), J, inj)


_Q.A.grid = _grid_exempt
_P4g = _Q.P4.grid


def _p4grid_own(J, K=2, *a, **k):
    """odd rows: the same own-row node rule as _drop_own (an injected own row never falls back to the K2 node at a dropped x)."""
    inj = {float(D) for D in V.INJ.get(J, [])}
    g = _P4g(J, K, *a, **k)
    return _drop_own(g, J, inj) if inj else g


_Q.P4.grid = _p4grid_own
