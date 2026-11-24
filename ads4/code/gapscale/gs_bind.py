"""
gapscale/gs_bind.py (GAP_SCALING.md §4): the binding state (J = 5, +-, parity-odd) at Delta_g = 2000 on its own-mass K2 row
(Delta = 2005, gs_rows.py), with the NEARGAP2 error budget run UNCHANGED (neargap2/ng_bind.run_state, read-only import), and the
H_mono test: Gamma(x) in c2 units of the J = 5 rows at Delta = 200, 305, 400, 800, 1005, 1600, 2005.
Two redirections only (in-process, nothing in neargap2/ or spinfunc/ is written; ng_guard, imported by ng_chain, would refuse it):
  * the row loader w2_chain._pk also reads gapscale/cache/rows (only when the chain is not in its NOOWN 'published rows' mode);
  * ng_bind.gauge_factor at Delta_g = 2000 (no u3/x2 pickles exist there) uses the min-norm weight of gs_weights.py
    (cache/weights_2000.npz): max(1, |Re X/psi2|) as a local max over +-0.5 on nu <= 25, as in ng_bind.
Output goes to stdout only (redirect it; ng_guard forbids writes here).
Usage: python3 gs_bind.py > cache/bind.txt
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
GROWS = os.path.join(HERE, 'cache', 'rows')
WTS = os.path.join(HERE, 'cache', 'weights_2000.npz')
sys.path.insert(0, os.path.join(HERE, '..', 'neargap2'))
import numpy as np
_w2000 = dict(np.load(WTS))                     # load before the guard matters (reads are allowed anyway)
import ng_chain as NC
import ng_bind as NB
W, Q, E = NC.W, NC.Q, NC.E
_ngpk = W._pk


def _pk(pats):
    R = _ngpk(pats)
    if not NC.NOOWN:
        for p in pats:
            for fn in sorted(glob.glob(os.path.join(GROWS, p))):
                R.update(pickle.load(open(fn, 'rb')))
    return R


W._pk = _pk
_gf_orig = NB.gauge_factor


def gauge_factor(Dg, nu):
    if int(round(Dg)) != 2000:
        return _gf_orig(Dg, nu)
    g = _w2000['nu']; X = _w2000['X'].real; phi = _w2000['phi']
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


def pr(s=''):
    print(s, flush=True)


def hmono():
    pr("## H_mono at J = 5 (+-, parity-odd): Gamma(x) in c2 = m^2 x relative correction, K = 2, on every J = 5 row")
    NC.NOOWN = False; NC.KK, NC.PK = 2, None
    g = Q.odd_grid(5)
    Ds = sorted(g)
    pr("   rows: " + ", ".join(f"Delta={D:g} (m={g[D]['m']:.1f}, J^2/m={25/g[D]['m']:.4f})" for D in Ds))
    xs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    for kd in ('C2', 'C3H'):
        pr(f"   {kd}:  x  " + "  ".join(f"{D:>8g}" for D in Ds))
        for x in xs:
            vals = []
            for D in Ds:
                xr = g[D]['xs']
                vals.append(f"{np.interp(x, xr, g[D][kd]) if xr[0] <= x <= xr[-1] else np.nan:+8.1f}")
            pr(f"        {x:.1f} " + "  ".join(vals))
    return g


if __name__ == '__main__':
    hmono()
    pr("\n## ng_bind.run_state (unchanged) at the binding state")
    for Dg, Down, mfs in ((1000, 1005.0, (1.0049925,)), (2000, 2005.0, (float(E.mgap(2000, 5)), 1.0049925, 1.01, 1.02, 1.05))):
        pr(f"\n# Delta_g = {Dg}: own row Delta = {Down:g}; twist edge m/Lam = {E.mgap(Dg, 5):.7f}")
        for mf in mfs:
            NB.run_state(Dg, 5, Down, mf, pr)
