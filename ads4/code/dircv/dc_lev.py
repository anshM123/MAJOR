"""
dc_lev (DIRCV.md §5): final estimates at (1100, 100) on a record set TAG (default the pooled '_pool'): the lambda family
(record x (exact/MC scalar)^lambda; 0 = dir, 1 = cv), the unbiased P_J-cv twin, the k = 3 check, lambda constraints and leverages.
Usage: DC_LOCAL=_pool,_poolP python3 dc_lev.py [_pool]  -> cache/dc_lev<TAG>.txt
"""
import dc_guard, os, sys
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
TAG = sys.argv[1] if len(sys.argv) > 1 else '_pool'
PTAG = TAG + 'P'
N4 = [13.3, 16.3, 20.8, 25.0]
N6 = [13.3, 16.3, 16.7, 20.8, 25.0, 33.4]
TB = {'++': (1158.0, 1441.0), '+-': (1152.0, 1439.0)}      # T_B with |Y| / |X| (N2_SUBREGGE §2.1, c = 3)
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)

def cfg(lam, tag, kf, est='cv'):
    LAM['v'] = lam
    return DP.run(J, D, est, [tag], kf, pr=lambda s: None)

pr(f"# dc_lev (1100,100) record set {TAG} (P_J twin {PTAG})")
Pcv = {kf: cfg(1.0, PTAG, kf) for kf in (False, True)}
Pdir = {kf: cfg(0.0, PTAG, kf) for kf in (False, True)}
models = (('F1', False, N4), ('F1', True, N6), ('F1B', True, N6))
res = {}
for lam in (0.0, 0.25, 0.5, 0.75, 1.0):
    for kf in (False, True):
        res[(lam, kf)] = cfg(lam, TAG, kf)
pr("## r (+- jackknife), rms, and paired r - r(P_J-cv twin)")
for model, kf, nus in models:
    for h in ('++', '+-'):
        fp = F.fit(Pcv[kf], h, nus, model); fpd = F.fit(Pdir[kf], h, nus, model)
        s = [f"PJcv {fp['r']:+6.0f}+-{fp['er']:4.0f} (rms {fp['rms']:.0f}); PJdir {fpd['r']:+6.0f}+-{fpd['er']:4.0f}"]
        for lam in (0.0, 0.5, 1.0):
            p = F.paired(res[(lam, kf)], Pcv[kf], h, nus, model)
            s.append(f"lam {lam}: {p['fa']['r']:+6.0f}+-{p['fa']['er']:4.0f} rms {p['fa']['rms']:3.0f}, -PJcv {p['dr']:+6.0f}+-{p['edr']:4.0f}")
        pr(f"  {model} {len(nus)} nodes {h}: " + " | ".join(s))
# lambda constraints
pr("## lambda constraints (lambda = graviton wave bias / scalar wave bias)")
cons = []
for model, kf, nus in models:
    for h in ('++', '+-'):
        d = [F.paired(res[(l, kf)], Pcv[kf], h, nus, model) for l in (0.0, 1.0)]
        # dr(lam) ~ dr0 + lam (dr1 - dr0); zero at lam* ; error from the jackknife of dr at lam*  (iterate once)
        lam0 = -d[0]['dr']/(d[1]['dr'] - d[0]['dr'])
        e = 0.5*(d[0]['edr'] + d[1]['edr'])/abs(d[1]['dr'] - d[0]['dr'])
        pr(f"  P_J-cv twin, {model} {len(nus)} nodes {h}: lambda = {lam0:+.2f} +- {e:.2f}")
        cons.append((model, len(nus), h, lam0, e))
for kf, nus in ((False, N4), (True, N6)):
    h, sg = '+-', -1.0
    vals = []
    for l in (0.0, 1.0):
        o = res[(l, kf)]; f2 = F.fit(o, h, nus, 'F1'); c3, c3j = c3fit(o, h, nus); R = c3j.shape[0]
        dj = c3j[:, 0] - sg*f2['cj'][:, 0]
        vals.append((c3[0] - sg*f2['c'][0], np.sqrt((R - 1)/R*np.sum((dj - dj.mean())**2))))
    lam0 = -vals[0][0]/(vals[1][0] - vals[0][0]); e = 0.5*(vals[0][1] + vals[1][1])/abs(vals[1][0] - vals[0][0])
    pr(f"  k=3 (+-), {len(nus)} nodes: diff(lam=0) {vals[0][0]:+.0f}+-{vals[0][1]:.0f}, diff(lam=1) {vals[1][0]:+.0f}+-{vals[1][1]:.0f} -> lambda = {lam0:+.2f} +- {e:.2f}")
# leverage table
pr("## leverage L = T_B/(|r| + rms) [|Y| .. |X|];  2 sigma: T_B^|Y|/(|r| + 2 sig + rms)")
for model, kf, nus in models:
    for name, o in [('dir (lam 0)', res[(0.0, kf)]), ('lam 0.5', res[(0.5, kf)]), ('cv (lam 1)', res[(1.0, kf)]), ('P_J-cv', Pcv[kf])]:
        s = []
        for h in ('++', '+-'):
            f = F.fit(o, h, nus, model); a = abs(f['r']) + f['rms']
            s.append(f"{h}: r {f['r']:+5.0f}+-{f['er']:3.0f} rms {f['rms']:3.0f}: L {TB[h][0]/a:.1f}-{TB[h][1]/a:.1f}, 2sig {TB[h][0]/(a + 2*f['er']):.1f}")
        pr(f"  {model} {len(nus)}n {name:12s} " + " | ".join(s))
open(os.path.join(dc_guard.CACHE, 'dc_lev%s.txt' % TAG), 'w').write('\n'.join(lines) + '\n')
