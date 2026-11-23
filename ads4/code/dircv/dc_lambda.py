"""
dc_lambda (DIRCV.md §3.6): one-parameter family of per-wave estimators, record x (exact/MC scalar)^lambda:  lambda = 0 is dir,
lambda = 1 is cv; lambda = graviton relative wave bias / scalar relative wave bias makes the estimator unbiased.
For each lambda: r (F1, 4 pw nodes; F1B, 6 nodes), paired difference to the unbiased P_J-cv twin (same CRN), k = 3 check (+-).
Usage: python3 dc_lambda.py [tag] -> cache/dc_lambda<tag>.txt
"""
import dc_guard, os, sys, pickle
import numpy as np
import z3_nodes as Z3
import dc_pipe as DP, dc_fit as F
from dc_k3 import c3fit
LAM = {'v': 1.0}
_sc = Z3.scalar_check


def scalar_check(*a, **k):
    o = dict(_sc(*a, **k))
    o['cv'] = np.exp(LAM['v']*np.log(o['cv'])); o['cv_jk'] = np.exp(LAM['v']*np.log(o['cv_jk']))
    return o


Z3.scalar_check = scalar_check
J, D = 100, 1100.0
tag = sys.argv[1] if len(sys.argv) > 1 else '_w'
N4 = [13.3, 16.3, 20.8, 25.0]
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)
LAM['v'] = 1.0
ref = DP.run(J, D, 'cv', [tag + 'P'] if tag != '_w' else ['_wP'], False, pr=lambda s: None)
pr(f"# lambda family, tag {tag}; P_J-cv twin r (F1, 4 nodes): " + ", ".join(f"{h}: {F.fit(ref, h, N4)['r']:+.0f} +- {F.fit(ref, h, N4)['er']:.0f}" for h in ('++', '+-')))
for lam in (0.0, 0.25, 0.5, 0.75, 1.0):
    LAM['v'] = lam
    o = DP.run(J, D, 'cv', [tag], False, pr=lambda s: None)
    s = []
    for h, sg in (('++', 1.0), ('+-', -1.0)):
        p = F.paired(o, ref, h, N4)
        f2 = p['fa']; c3, c3j = c3fit(o, h, N4); R = c3j.shape[0]
        dj = c3j[:, 0] - sg*f2['cj'][:, 0]; e = np.sqrt((R - 1)/R*np.sum((dj - dj.mean())**2))
        s.append(f"{h}: r = {f2['r']:+6.0f} +- {f2['er']:4.0f}; r - r(PJcv) = {p['dr']:+6.0f} +- {p['edr']:4.0f}; k3 diff {c3[0] - sg*f2['c'][0]:+6.0f} +- {e:5.0f}")
    pr(f"  lambda={lam:.2f}: " + " | ".join(s))
open(os.path.join(dc_guard.CACHE, 'dc_lambda%s.txt' % tag), 'w').write('\n'.join(lines) + '\n')
