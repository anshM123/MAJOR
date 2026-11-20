"""
ee1000_shapecheck (EE_DG1000_LARGEJ.md §4): validation of the shape-fit representation of Gamma(x) on a row: (a) residuals, (b) cross-
validation (fit on every other node, predict the others), (c) the gap-state dF per x-bin for spline vs shape (K = 3..6) and for the two
half-node fits, (d) the Gamma table.  Usage: W2_ROWS=<Delta> python3 ee1000_shapecheck.py J  -> cache/ee1000_shapecheck_J<J>_<Delta>.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import ee1000_chain as CH
E, Q = CH.E, CH.Q
J = int(sys.argv[1])
odd = J % 2 == 1
g = CH.trim(Q.odd_grid(J) if odd else Q.even_grid(J))
D = sorted(g)[0]; r = g[D]; xs = r['xs']; m = r['m']
keys = ['C2', 'C3H'] if odd else [(kd, h) for kd in ('C2', 'C3H') for h in ('++', '+-')]
hs = ('+-',) if odd else ('++', '+-')
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)
lab = lambda k: k if isinstance(k, str) else k[0] + ' ' + k[1]
pr(f"# J={J}, Delta={D:.0f}, m={m:.1f}, {len(xs)} nodes (after Hampel), x = {xs[0]:.4f} ... {xs[-1]:.4f}; despiked: {CH.W.DESPIKE_LOG}")
for k in keys:
    pr(f"  Gamma {lab(k):7s}: " + " ".join(f"{x:.4f}:{v:.0f}" for x, v in zip(xs, r[k])))
def subrow(idx):
    rr = dict(r); rr['xs'] = xs[idx]
    for k in keys:
        rr[k] = np.asarray(r[k])[idx]
    return rr
for K in (3, 4, 5, 6):
    sf = CH.ShapeFit(J, r, keys, K=K)
    ev, od = np.arange(0, len(xs), 2), np.arange(1, len(xs), 2)
    cv = []
    for k in keys:
        e1 = CH.ShapeFit(J, subrow(ev), [k], K=K)(k, xs[od]) - np.asarray(r[k])[od]
        e2 = CH.ShapeFit(J, subrow(od), [k], K=K)(k, xs[ev]) - np.asarray(r[k])[ev]
        cv.append(f"{lab(k)}: res {sf.res[k][0]/sf.res[k][1]:.3f}, CV {np.sqrt(np.mean(np.r_[e1, e2]**2))/sf.res[k][1]:.3f}")
    pr(f"  K={K}: " + "; ".join(cv))
bins = [0, xs[0], 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, xs[-1], 1.0]
mfl = E.mgap(1000, J)
st = {h: E.hyb(1000, J, mfl, h, 256) for h in hs}
def contrib(fun):
    out = {}
    for h in hs:
        s = st[h]; x = s['x']
        k2, k3 = ('C2', 'C3H') if odd else (('C2', h), ('C3H', h))
        env = Q.P4.env if odd else Q.A.env
        c = (s['p2']*env(h, x, 2)*fun(k2, x) + s['p3']*s['u']*env(h, x, 3)*fun(k3, x))/m**2/s['Ns']
        out[h] = ([c[(x >= bins[i]) & (x < bins[i + 1])].sum() for i in range(len(bins) - 1)], c.sum())
    return out
from scipy.interpolate import CubicSpline
spl = {k: CubicSpline(xs, r[k], bc_type='natural') for k in keys}
variants = [('spline', lambda k, x: spl[k](np.minimum(x, xs[-1])))]
for K in (3, 4, 5, 6):
    sf = CH.ShapeFit(J, r, keys, K=K); variants.append((f'shape K={K}', sf))
for nm, idx in (('even nodes', np.arange(0, len(xs), 2)), ('odd nodes', np.arange(1, len(xs), 2))):
    sf = CH.ShapeFit(J, subrow(idx), keys, K=4); variants.append((f'shape K=4 {nm}', sf))
pr(f"# gap state m/Lam = {mfl:.4f}: dF_hh/N_h per x-bin " + " ".join(f"[{bins[i]:.3g},{bins[i+1]:.3g})" for i in range(len(bins) - 1)))
for nm, fun in variants:
    c = contrib(fun)
    for h in hs:
        pr(f"  {nm:22s} {h}: " + " ".join(f"{v:+.3f}" for v in c[h][0]) + f"   total {c[h][1]:+.3f}")
open(os.path.join(CH.C, 'ee1000_shapecheck_J%d_%d.txt' % (J, D)), 'w').write("\n".join(lines) + "\n")
