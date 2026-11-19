"""
r2_ladder (LARGEJ_OFFDIAG.md §2): the real helicity off-diagonal Gamma_off along LADDERS OF FIXED J/m (fixed AdS impact parameter beta,
i.e. fixed b-hat = (2J+1) Lam/m at fixed Delta_g: the leading-Regge direction m -> oo at fixed (nu, beta)), compared at equal nu.
Uses q3_analysis.even_grid unchanged (C^_2', fixed C^_3, phase strip, row alignment) with the r1 caches (r1_run.py) merged into the q2
loader, and optionally the direct (RQMC) rows of q4 / r1 (env R2_DIRECT=1).
Quantities (l9 units): F_+- density per unit psi-weight = eo2 * Gamma_off / m  (q3 / l9_neargap), so at fixed nu the off-diagonal
density scales as Gamma_off/m.  Scalings tested at fixed (nu, J/m):  Gamma_off ~ m^s  ->  density ~ m^(s-1).
   s = 1 : m-independent density  = a LEADING-REGGE off-diagonal (would contradict helicity-diagonal leading Regge)
   s < 1 : off-diagonal vanishes in the Regge limit; s = -1 is "O(1/m^2) at fixed (nu, beta)".
Usage: python3 r2_ladder.py  -> cache/r2_ladder.txt
"""
import os, sys, glob, pickle, io
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.interpolate import CubicSpline
import q3_analysis as Q

CACHE = os.path.join(HERE, 'cache')
_orig_q2 = Q._q2_even


def _q2_even(kind, J):
    R = _orig_q2(kind, J)
    for fn in sorted(glob.glob(os.path.join(CACHE, 'r1_pts_%s_J%d_k2.pkl' % (kind, J)))):
        R.update(pickle.load(open(fn, 'rb')))
    return R


Q._q2_even = _q2_even


def rows(J):
    g = Q.even_grid(J)
    return g


LADDERS = [
    [(6, 400.0), (12, 800.0), (24, 1600.0)],
    [(8, 400.0), (16, 800.0), (32, 1600.0)],
    [(10, 400.0), (20, 800.0)],
    [(14, 800.0), (28, 1600.0)],
    [(18, 336.0), (36, 672.0), (72, 1344.0)],
]


def main():
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    cache = {}
    def get(J):
        if J not in cache:
            try:
                cache[J] = rows(J)
            except Exception as ex:
                cache[J] = {}; pr(f"  J={J}: no grid ({ex})")
        return cache[J]
    pr("# r2: Gamma_off (C2, C3H; aligned, l9 units) along fixed-J/m ladders, at common nu (cubic spline in nu on each row)")
    summary = []
    for lad in LADDERS:
        have = [(J, D) for J, D in lad if D in get(J)]
        if len(have) < 2:
            continue
        recs = [get(J)[D] for J, D in have]
        ms = [r['m'] for r in recs]
        pr(f"\n## ladder J/m = " + ", ".join(f"({J},{D:.0f}): {J/m:.4f}, m={m:.0f}, J^2/m={J*J/m:.2f}" for (J, D), m in zip(have, ms)))
        numax = min(r['xs'][-1]*r['m'] for r in recs); numin = max(max(r['xs'][0]*r['m'], 15.0) for r in recs)
        nus = np.linspace(numin, min(numax, 400.0), 25)
        for kd in ('C2', 'C3H'):
            spl = [CubicSpline(r['xs']*r['m'], r[(kd, 'off')]) for r in recs]
            vals = np.array([s(nus) for s in spl])
            # common sign per row: maximise overlap with the first row
            for i in range(1, len(vals)):
                if vals[i] @ vals[0] < 0:
                    vals[i] = -vals[i]
            pr(f"  {kd} Gamma_off at nu:  " + " ".join(f"{v:6.0f}" for v in nus))
            for (J, D), v in zip(have, vals):
                pr(f"     ({J:2d},{D:5.0f})          " + " ".join(f"{a:6.2f}" for a in v))
            # scaling exponent s from the L2 norm over the common nu window and the peak
            nrm = np.sqrt((vals**2).mean(1)); pk = np.abs(vals).max(1)
            for i in range(1, len(vals)):
                s_n = np.log(nrm[i]/nrm[0])/np.log(ms[i]/ms[0]); s_p = np.log(pk[i]/pk[0])/np.log(ms[i]/ms[0])
                pr(f"     scaling ({have[0][0]},{have[0][1]:.0f}) -> ({have[i][0]},{have[i][1]:.0f}): Gamma_off ~ m^s, s = {s_n:+.2f} (rms), "
                   f"{s_p:+.2f} (peak)  => off density ~ m^{s_n - 1:+.2f}")
                summary.append((kd, have[0], have[i], s_n, s_p))
            # also the diagonal corrections for reference (C2 ++), same window
        for kd in ('C2',):
            for h in ('++', '+-'):
                spl = [CubicSpline(r['xs']*r['m'], r[(kd, h)]) for r in recs]
                vals = np.array([s(nus) for s in spl])
                nrm = np.sqrt((vals**2).mean(1))
                pr(f"  {kd} {h} Gamma_hh rms on the window: " + ", ".join(f"({J},{D:.0f}) {v:.0f}" for (J, D), v in zip(have, nrm)) +
                   "   (density ~ Gamma_hh/m^2)")
    pr("\n# summary (kind, from, to, s_rms, s_peak)")
    for s in summary:
        pr(f"  {s[0]:3s} {s[1]} -> {s[2]}: s = {s[3]:+.2f} / {s[4]:+.2f}")
    open(os.path.join(CACHE, 'r2_ladder%s.txt' % os.environ.get('R2_TAG', '')), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
