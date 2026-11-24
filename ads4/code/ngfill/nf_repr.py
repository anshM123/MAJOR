"""
nf_repr (NEARGAP_FILL.md §3): a noise-aware representation term for the odd-J direct-bulk twist edges (J = 23, 31; J = 27 as control).
The protocol (highspin/hs_odd.py + smallnu/sn_odd.py, both run UNCHANGED via runpy / import, on extfold's merged point set 'all',
read-only) fits the direct bulk nodes Gamma(x_i) with the shape basis c(x)D_h + s(x)xD_h'/(2J+1) + r(x) (degree K = 6) by ORDINARY
least squares, and sets E_repr = 3 max|lam(K6) - lam(K4, K5, K7, two half-samples at K4)|.  On a direct row the node errors are
heteroscedastic (sigma_i grows ~10x from x = 0.05 to x = 0.8), so OLS with 3(K+1) = 24 parameters (K7) on 34 nodes lets the noisy
large-x nodes drive the fit: the K7 probe is then a noise amplifier (EXT_FOLD §3.2: its jackknife sigma is 5x that of K6).
Mode 'wls' replaces ONLY the estimator: every shape fit (main fit, all probes, all jackknife replicas, sn_odd's bulk model) becomes
weighted least squares with weights 1/sigma_i^2, sigma_i^2 = mean over the two kinds of the node's replica-jackknife variance
(computed once, from the full point set; held fixed in the jackknife).  The rule E_repr = 3 max|...| over the SAME probes is kept.
Mode 'ols' reruns the protocol unchanged (reproduction of EXT_FOLD's numbers).
Diagnostics (both modes): per degree K = 4..8, lam and its jackknife sigma, and the fit's chi^2/dof (per kind, node sigma_i).
Usage: python3 nf_repr.py J mode [mfl,...]   -> cache/repr_<mode>/hs_odd_J<J>.{txt,json}, sn_odd_J<J>.{txt,json}, nf_repr_J<J>.txt
"""
import os, sys, json, runpy
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import nf_guard  # noqa: F401  (first; stands in for ef_guard)
sys.path.insert(0, os.path.join(HERE, '..', 'extfold'))
import ef_core as C              # read-only (merged point sets of EXT_FOLD)
import numpy as np
H = C.H
OLS_FIT = H.ShapeFit
SIG = {}


class WShapeFit(OLS_FIT):
    """hs_core.ShapeFit with weights 1/sigma_i^2 (SIG, keyed by round(x, 6)); falls back to OLS if any node has no sigma."""
    def __init__(self, J, xs, y, h, K=4, x0=0.0):
        xs = np.asarray(xs, float); y = np.asarray(y, float)
        s = np.array([SIG.get(round(float(x), 6), np.nan) for x in xs])
        if not np.all(np.isfinite(s)) or x0 > 0:
            return OLS_FIT.__init__(self, J, xs, y, h, K=K, x0=x0)
        self.J, self.K, self.h, self.x0 = J, K, h, x0
        self.xm = float(np.max(xs))
        A = H._basis(h, J, xs, self.xm, K, x0)
        Aw = A/s[:, None]; sc = np.abs(Aw).max(0); sc[sc == 0] = 1
        c, *_ = np.linalg.lstsq(Aw/sc, y/s, rcond=None)
        self.c = c/sc
        res = y - A @ self.c
        self.rms_rel = float(np.sqrt(np.mean(res**2))/max(np.sqrt(np.mean(y**2)), 1e-30))
        self.maxres = float(np.abs(res).max())
        self.chi2 = float(((res/s)**2).sum()); self.dof = len(y) - 3*(K + 1)


if __name__ == '__main__':
    J, mode = int(sys.argv[1]), sys.argv[2]
    D = 1000.0 + J
    C.SET[0] = 'all'
    out = os.path.join(HERE, 'cache', 'repr_' + mode); os.makedirs(out, exist_ok=True)
    H.HSC = out; C.S.SNC = out
    import hs_odd as HO          # read-only import
    # node sigmas: replica jackknife of every direct node, both kinds
    reps = [HO.rows(J, D, jk=j)[1] for j in range(8)]
    rK, rI, isdir = HO.rows(J, D)
    sk = {}
    for i, x in enumerate(rI['xs']):
        if not isdir[i]:
            continue
        v = []
        for kd in ('C2', 'C3H'):
            vj = np.array([rr[kd][int(np.argmin(np.abs(rr['xs'] - x)))] for rr in reps])
            v.append(HO.jk_sigma(vj)); sk[(kd, round(float(x), 6))] = v[-1]
        SIG[round(float(x), 6)] = float(np.sqrt(np.mean(np.square(v))))
    if mode == 'wls':
        H.ShapeFit = WShapeFit
    lines = []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    mfls = sys.argv[3] if len(sys.argv) > 3 else None
    sys.argv = ['hs_odd.py', str(J)] + ([mfls] if mfls else [])
    runpy.run_path(os.path.join(HERE, '..', 'highspin', 'hs_odd.py'), run_name='__main__')
    C.H.HSC = out
    import sn_odd as SO          # read-only import; reads out/hs_odd_J<J>.json (FULL normalisation, as EXT_FOLD)
    SO.main(J)
    # diagnostics: degree scan with jackknife, chi^2/dof per kind (node sigmas as above)
    bad = H.refbad(J, D, HO.RTOL)
    rowm = H.HN.mphys(D, J); mtw = rowm/H.E.lam_of(1000.0); HO.XROW[0] = 1/mtw
    dset = set(np.round(H.direct_x(J, D), 6)) - H.failed_x(J, D)
    isd = lambda rr: np.array([round(x, 6) in dset for x in rr['xs']])
    sel = (rI['xs'] >= 0.05 - 1e-9) & (rI['xs'] < HO.XCUT - 1e-9) & isdir
    pr(f"# nf_repr J={J} mode={mode}: {int(sel.sum())} direct bulk nodes; node sigma (c2) median {np.median([SIG[round(float(x), 6)] for x in rI['xs'][sel]]):.0f}, "
       f"x<0.3 {np.median([SIG[round(float(x), 6)] for x in rI['xs'][sel] if x < 0.3]):.0f}, x>0.6 {np.median([SIG[round(float(x), 6)] for x in rI['xs'][sel] if x > 0.6]):.0f}")
    xb = rI['xs'][sel]
    for K in (4, 5, 6, 7, 8):
        s2 = []
        for kd in ('C2', 'C3H'):
            yb = np.asarray(rI[kd], float)[sel]; sb = np.array([sk[(kd, round(float(x), 6))] for x in xb])
            f = H.ShapeFit(J, xb, yb, '+-', K=K)
            r = yb - f(xb); s2.append(f"{kd} chi2/dof {((r/sb)**2).sum():.1f}/{len(yb) - 3*(K + 1)}")
        pr(f"   K={K}: " + "; ".join(s2))
    res = {}
    for mfl in [mtw] + ([float(v) for v in mfls.split(',')] if mfls else []):
        o = {}
        for K in (4, 5, 6, 7, 8):
            M = HO.Model(J, rK, rI, isdir, bad, bulk='dir', kb=K)
            o[K] = HO.value(J, mfl, M, 'gk')[0]
            lj = [HO.value(J, mfl, HO.Model(J, rK, rr, isd(rr), bad, bulk='dir', kb=K), 'gk')[0] for rr in reps]
            dj = [l - HO.value(J, mfl, HO.Model(J, rK, rr, isd(rr), bad, bulk='dir', kb=6), 'gk')[0] for l, rr in zip(lj, reps)]
            o[f'jk{K}'] = HO.jk_sigma(lj); o[f'djk{K}'] = HO.jk_sigma(dj)
        res[f'{mfl:.7f}'] = o
        pr(f"  m/Lam={mfl:.7f}: " + ", ".join(f"K{K} {o[K]:+.4f} (jk {o[f'jk{K}']:.4f}; K-K6 {o[K] - o[6]:+.4f} +- {o[f'djk{K}']:.4f})" for K in (4, 5, 6, 7, 8)))
    open(os.path.join(out, 'nf_repr_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(out, 'nf_repr_J%d.json' % J), 'w'), indent=1, default=float)
