"""
Independent validation of the LC spinning tensors for l = 2 (and l = 1): evaluate, with setting A, all radial-series-converged e2e grid points
with u < 0.25 and compare with the cached radial values; report the deviation as a function of t and its weighted impact on the moments,
sum_p |LC_p - rad_p| w_p / |moment| (the error the LC representation would cause if used at those points), for the modes (l,-l), (l,l).
Usage: python3 v2_l2_validation.py       Output: out_v2_l2_validation.txt
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np
from moments import TensorMoments
from corner import evaluate_points, e2e_file, PAR_A, TensorMomentsLC
from t4_spinning_e2e import main_cases
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
for (l, J) in ((2, 8), (2, 20), (1, 20)):
    Delta = [c[0] for c in main_cases() if c[1] == J][0]
    base = TensorMoments.load(e2e_file(l, Delta, J)); g = base.g
    sel = g['ok'] & base.conv & (g['U'] < 0.25)
    (ls, fm, it, info), = evaluate_points(l, Delta, J, 1.75, g['U'][sel], g['S'][sel], pars=(PAR_A,))
    T = 2*g['S']*(Delta - J - 2 + (Delta + J - 1)*g['U'])
    pr(f"# l={l} (Delta,J)=({Delta:.3f},{J}): {int(sel.sum())} converged points with u < 0.25")
    keys = sorted(base.modes.keys())
    A = np.array([base.modes[k][sel]*np.exp(base.logscale[sel] - ls) for k in keys]); B = np.array([fm[k] for k in keys])
    dev = np.max(np.abs(A - B), axis=0)/np.max(np.abs(A), axis=0)
    for (t0, t1) in ((0, 2), (2, 4), (4, 8), (8, 16), (16, 90)):
        m = (T[sel] >= t0) & (T[sel] < t1)
        if m.any():
            pr(f"   t in [{t0},{t1}): {int(m.sum())} pts, max dev {dev[m].max():.1e}, median {np.median(dev[m]):.1e}")
    # weighted impact on the moments
    mom = TensorMomentsLC.__new__(TensorMomentsLC)
    mom.l, mom.Delta, mom.J, mom.Dphi = l, Delta, J, 1.75
    mom.g = g; mom.logscale = base.logscale; mom.modes = base.modes
    dm = {k: np.zeros_like(base.modes[k]) for k in keys}
    for j, k in enumerate(keys):
        dm[k][sel] = np.abs(A[j] - B[j])*np.exp(ls - base.logscale[sel])
    for k in (2, 3, 4):
        row = []
        for q in ((l, -l), (l, l)):
            dens, sh = mom._density(k, 0.0, *q); de, _ = mom._density(k, 0.0, *q, modes=dm)
            row.append(f"{q}: {np.abs(de).sum()/abs(dens.sum()):.1e}")
        pr(f"   k={k}, nu=0: weighted |LC - radial| / |moment| over these points: " + "  ".join(row))
open(os.path.join(HERE, 'out_v2_l2_validation.txt'), 'w').write("\n".join(out) + "\n")
