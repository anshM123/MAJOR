"""
sn_core (SMALLNU_WINDOW.md): Gamma at small nu (nu < 50) on a near-gap state's row, in EXACTLY the chain's convention, from direct
points (sn_points.py, smallnu/cache) -- through highspin/hs_core (imported read-only; its rows = the chain's K2 rows with direct
nodes injected through neargap2/ng_val's loader wrappers, hc_nodes records, kappa image for odd J).

* Direct points are read from highspin/cache (HIGHSPIN's rows) AND smallnu/cache (this note; smallnu wins on a common key):
  hc_nodes.pts is replaced by a merging reader.  BASE[0] = True hides the smallnu points (the 'baseline' grid = HIGHSPIN's).
* Normalisation.  Each row's per-state constant a (odd: a_{+-}; even: a_{++}, a_{+-}) is fitted on the row's C2 nodes x <= 0.25,
  nu >= 15 (p4_odd_analysis.grid / l9_analysis.grid), and every record carries exp(-c0), c0 = the log-ratio at the row's first node.
  New small-nu nodes enter that fit (and become the first node), which moves the effective normalisation by ~1e-4 and so shifts Gamma
  on the WHOLE row by ~1e-4 m^2 S(x) ~ 1e2 c2 (a near-pure renormalisation of the image).  To report the new nodes in the SAME
  normalisation as the baseline (HIGHSPIN) row, the injected grid is transformed back (_transform): one ratio r_h per helicity measured
  on the old nodes, Gamma' = m^2 [(Gamma/m^2 + S)/r_h - S], S = sh_k/env_k + sg_k D_h (the chain module's own functions);
  off-diagonal / sqrt(r_++ r_+-).  Checked on the old nodes (x < 0.94; reproduces the baseline to rounding).
* Jackknife: replica jk left out in BOTH grids (all injected nodes together), each transformed with its own baseline a.
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import sn_guard  # noqa: F401  (first; also stands in for hs_guard / ng_guard)
HS = os.path.join(HERE, '..', 'highspin')
sys.path.insert(0, HS); sys.path.insert(0, HERE)
import numpy as np
import hs_core as H
HN, Q, E = H.HN, H.Q, H.E
SNC = os.path.join(HERE, 'cache')
HSC = H.HSC
BASE = [False]
SUBDIR = [None]          # read the smallnu points from smallnu/cache/<subdir> instead (cross-check caches)


def _pts(kind, J, D):
    f = 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)
    out = {}
    p = os.path.join(HSC, f)
    if os.path.exists(p):
        out.update(pickle.load(open(p, 'rb')))
    if not BASE[0]:
        q = os.path.join(SNC if SUBDIR[0] is None else os.path.join(SNC, SUBDIR[0]), f)
        if os.path.exists(q):
            out.update(pickle.load(open(q, 'rb')))
    return out


HN.pts = _pts
HN.pfile = lambda kind, J, D: os.path.join(HSC, 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D))


def sn_nus(J, D, sub=None):
    """nu of the smallnu direct points of state D (state AND reference present)."""
    f = os.path.join(SNC if sub is None else os.path.join(SNC, sub), 'hc_pts_EE_J%d_D%g.pkl' % (J, D))
    if not os.path.exists(f):
        return []
    P = pickle.load(open(f, 'rb'))
    return sorted(k[1] for k in P if k[0] == D and (H.DREF, k[1]) in P)


def rel_err(J, D, nu, sub=None):
    """max relative replica error of the state and reference records (both kinds) at nu."""
    e = 0.0
    for kind in ('EE', 'C3H'):
        f = os.path.join(SNC if sub is None else os.path.join(SNC, sub), 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D))
        if not os.path.exists(f):
            continue
        P = pickle.load(open(f, 'rb'))
        for DD in (D, H.DREF):
            for k, r in P.items():
                if k[0] == DD and abs(k[1] - nu) < 1e-4*max(nu, 1.0):
                    e = max(e, float(np.max(np.abs(r['err'])/np.maximum(np.abs(r['E']), 1e-300))))
    return e


def _rowgrid(J, D, base, jk=None, K=2):
    BASE[0] = base
    try:
        g = H.grid(J, K, inject=[D], jk=jk)
    finally:
        BASE[0] = False
    if J % 2:
        g = H.subtract_kappa(g, J)
    return g[D]


def _S(J, key, xs, m):
    """sh_k/env_k + sg_k D_h at xs for a diagonal key (the chain module's own functions)."""
    if J % 2:
        Mo, h, kd = Q.P4, '+-', key
        Dh = Mo.Dc
    else:
        Mo, (kd, h) = Q.A, key
        Dh = Mo.Dsh
    k = 2 if kd == 'C2' else 3
    sg = Mo.SIG[h] if k == 3 else 1.0
    return Mo.shift(J, h, k, xs, m)/Mo.env(h, xs, k) + sg*Dh(h, J, xs)


def _transform(J, rf, rb):
    """rf (row record of the injected grid) re-expressed in the normalisation of rb (baseline grid of the same row).
    Gamma/m^2 + S = g/a with g = G_raw/env (G_raw carries exp(-c0), c0 = the log-ratio at the row's FIRST node, which moves when a
    new node becomes the first one; it cancels in Gamma but not in a).  So the normalisation change is one number per helicity,
    r_h = (g/a)_f/(g/a)_b, measured on the common (old) nodes; new values: Gamma' = m^2 [(Gamma/m^2 + S)/r_h - S];
    off-diagonal: Gamma'_off = Gamma_off/sqrt(r_++ r_+-).  Returns (record, r, chk) with chk = max deviation (c2) on the old nodes."""
    out = dict(rf)
    m, xs = rf['m'], rf['xs']
    com = [(i, int(np.argmin(np.abs(rb['xs'] - x)))) for i, x in enumerate(xs) if np.min(np.abs(rb['xs'] - x)) < 1e-7]
    fi = np.array([c[0] for c in com]); bi = np.array([c[1] for c in com])
    hel = ['+-'] if J % 2 else ['++', '+-']
    r, chk = {}, 0.0
    for h in hel:
        Qs = {}
        for kd in ('C2', 'C3H'):
            key = kd if J % 2 else (kd, h)
            Sf = _S(J, key, xs, m); Sb = _S(J, key, rb['xs'], m)
            Qs[kd] = (np.asarray(rf[key])/m**2 + Sf, np.asarray(rb[key])/m**2 + Sb, Sf)
        rr = Qs['C2'][0][fi]/Qs['C2'][1][bi]
        r[h] = float(np.median(rr))
        for kd in ('C2', 'C3H'):
            key = kd if J % 2 else (kd, h)
            Qf, Qb, Sf = Qs[kd]
            out[key] = m*m*(Qf/r[h] - Sf)
            ok = xs[fi] < 0.94                    # (failed RQMC edge nodes, |Gamma| ~ 1e5: rounding only)
            chk = max(chk, float(np.max(np.abs(out[key][fi] - np.asarray(rb[key])[bi])[ok])))
    if not J % 2 and ('C2', 'off') in rf:
        f = np.sqrt(r['++']*r['+-'])
        for kd in ('C2', 'C3H'):
            out[(kd, 'off')] = np.asarray(rf[(kd, 'off')])/f
            chk = max(chk, float(np.max(np.abs(out[(kd, 'off')][fi] - np.asarray(rb[(kd, 'off')])[bi])[xs[fi] < 0.94])*m))
    return out, r, chk


def keys(J):
    return ['C2', 'C3H'] if J % 2 else [(kd, h) for kd in ('C2', 'C3H') for h in ('++', '+-', 'off')]


def small_nodes(J, D, numax=80.0, jkn=8):
    """Gamma at every node of row D with nu <= numax: per node dict(x, nu, src ('sn' = this note, 'hs' = HIGHSPIN direct, 'K2'),
    key -> (value, jackknife sigma)); values in the BASELINE normalisation.  Also returns the max transform residual on old nodes."""
    rb = _rowgrid(J, D, True); rf = _rowgrid(J, D, False)
    rt, rnorm, chk = _transform(J, rf, rb)
    reps = []
    for jk in range(jkn):
        b = _rowgrid(J, D, True, jk); f = _rowgrid(J, D, False, jk)
        t_, _, c_ = _transform(J, f, b)
        reps.append(t_); chk = max(chk, c_)
    m = rf['m']
    snx = {round(v/m, 6) for v in sn_nus(J, D)}
    hsx = set(np.round(H.direct_x(J, D), 6)) if True else set()
    out = []
    for i, x in enumerate(rt['xs']):
        if x*m > numax:
            continue
        rx = round(float(x), 6)
        src = 'sn' if rx in snx else ('hs' if rx in hsx else 'K2')
        rec = dict(x=float(x), nu=float(x*m), src=src)
        for k in keys(J):
            if k not in rt:
                continue
            v = float(rt[k][i])
            vj = []
            for rr in reps:
                j = int(np.argmin(np.abs(rr['xs'] - x)))
                vj.append(float(rr[k][j]) if abs(rr['xs'][j] - x) < 1e-7 else np.nan)
            vj = np.array(vj)
            e = float(np.sqrt((jkn - 1)/jkn*np.nansum((vj - np.nanmean(vj))**2))) if src != 'K2' else 0.0
            if k[-1] == 'off' if isinstance(k, tuple) else False:
                v, e = v*m, e*m                     # c_off = m Gamma_off (the quantity the chain integrates, OFFSCALE)
            rec[k] = (v, e)
        out.append(rec)
    da = {h: v - 1 for h, v in rnorm.items()}
    da = da['+-'] if J % 2 else da
    return out, chk, da, m
