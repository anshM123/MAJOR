"""
e5: the truncated CHMRS series (OMEGA_U / remk.py) against the exact moments on the e1 states.
Per skeleton and k: on series-'controlled' nodes, max |rho_ex - rho_ser| and how often it exceeds the series' own error estimate;
on 'uncontrolled' nodes, the size of rho_ex and of rho_ex - rho_ser; agreement in bands of x.
Usage: python3 e5_compare.py Dg  -> cache/e5_<Dg>.txt
"""
import os, sys, pickle, glob
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
Dg = int(sys.argv[1])
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)
for fn in sorted(glob.glob(os.path.join(HERE, 'cache', 'e1_%d_Df*.pkl' % Dg))):
    db = pickle.load(open(fn, 'rb'))
    pr(f"## {os.path.basename(fn)}: {len(db)} states")
    for k in range(2, 8):
        X, D, ERR, C, RE = [], [], [], [], []
        for r in db.values():
            rho_s, err, ctrl = r['ser'][k]
            X.append(r['x']); D.append(r['rho'][k] - rho_s); ERR.append(err); C.append(ctrl); RE.append(r['rho'][k])
        X, D, ERR, C, RE = map(np.concatenate, (X, D, ERR, C, RE))
        viol = C & (np.abs(D) > ERR)
        s = (f" k={k}: controlled nodes {C.sum()}: max|d rho| {np.abs(D[C]).max():.2e}, exceeds own error estimate on {viol.sum()} "
             f"(max excess ratio {np.max(np.abs(D[viol])/np.maximum(ERR[viol], 1e-300)) if viol.any() else 0:.3g}; at x up to {X[viol].max() if viol.any() else 0:.3f}); "
             f"uncontrolled {(~C).sum()}: max|rho_ex| {np.abs(RE[~C]).max() if (~C).any() else 0:.3g}, max|d rho| {np.abs(D[~C]).max() if (~C).any() else 0:.3g}")
        pr(s)
        v6 = viol & (X <= 0.6)
        pr(f"      violations at x <= 0.6: {v6.sum()} (max |d rho| there {np.abs(D[v6]).max() if v6.any() else 0:.2e}, max own estimate {ERR[v6].max() if v6.any() else 0:.2e})")
        row = "      max|d rho| by x band:"
        for lo, hi in ((0, .2), (.2, .3), (.3, .4), (.4, .5), (.5, .6), (.6, .7), (.7, .8), (.8, .9), (.9, 1.0)):
            sel = (X >= lo) & (X < hi)
            if sel.any():
                row += f" [{lo:.1f},{hi:.1f}) {np.abs(D[sel]).max():.1e}"
        pr(row)
open(os.path.join(HERE, 'cache', 'e5_%d.txt' % Dg), 'w').write("\n".join(out) + "\n")
