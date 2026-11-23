"""dc_nodes2 (DIRCV.md §8): per-node y' = Gamma - kappa - (A-1)kappa (+- jackknife) at (1100,100) for every estimator, all window nodes."""
import dc_guard, os, sys
import numpy as np
import z3_nodes as Z3
import dc_pipe as DP, dc_fit as F
LAM = {'v': 1.0}
_sc = Z3.scalar_check
def scalar_check(*a, **k):
    o = dict(_sc(*a, **k)); o['cv'] = np.exp(LAM['v']*np.log(o['cv'])); o['cv_jk'] = np.exp(LAM['v']*np.log(o['cv_jk'])); return o
Z3.scalar_check = scalar_check
TAG = sys.argv[1]
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)
def cfg(lam, tag, est='cv'):
    LAM['v'] = lam; return DP.run(100, 1100.0, est, [tag], True, pr=lambda s: None)
C = {'dir': cfg(0.0, TAG), 'lam0.7': cfg(0.7, TAG), 'cv': cfg(1.0, TAG), 'PJcv': cfg(1.0, TAG + 'P'), 'PJdir': cfg(1.0, TAG + 'P', 'dir')}
nus = [v for v in C['dir']['nus'] if 10.5 <= v <= 34.01]
pr(f"# per-node y' (C2), record set {TAG}; nodes {np.round(nus, 2).tolist()}; a_h: " + "; ".join(f"{k}: {o['a']}" for k, o in C.items()))
for h in ('++', '+-'):
    for k, o in C.items():
        f = F.fit(o, h, nus, 'F1')
        pr(f"  {h} {k:6s}: " + ", ".join(f"{v:+.0f}({e:.0f})" for v, e in zip(f['y'], f['ey'])) + f"   | F1 resid: " + ", ".join(f"{v:+.0f}" for v in f['res']))
open(os.path.join(dc_guard.CACHE, 'dc_nodes2%s.txt' % TAG), 'w').write('\n'.join(lines) + '\n')
