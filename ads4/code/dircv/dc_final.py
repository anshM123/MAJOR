"""
dc_final (DIRCV.md §0/§5): corrected remainder r and H_LJ(1000) leverage at (1100, 100) on a record set (default the pooled '_pool':
seeds 4242 + 1001 + 2002 at nu = 13.3, 16.3; seed 4242 elsewhere), with the per-wave estimator corrected by lambda = 0.5 +- 0.25
(graviton wave bias / scalar wave bias, §3), and the assumption-free P_J-cv twin as a cross-check.
Error budget: sigma_stat (jackknife at lambda = LAM0) (+) sigma_lambda = DLAM |r(1) - r(0)|.
Leverage: L = T_B/(|r| + rms) (T_B^|Y| .. T_B^|X|);  2 sigma: T_B^|Y|/(|r| + 2 sigma_tot + rms).  k = 3 check per estimator (+-).
Usage: DC_LOCAL=_pool,_poolP python3 dc_final.py [_pool] -> cache/dc_final<TAG>.txt
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
PTAG = '_wP' if TAG == '_w' else TAG + 'P'
LAM0, DLAM = 0.5, 0.25
N4 = [13.3, 16.3, 20.8, 25.0]; N6 = [13.3, 16.3, 16.7, 20.8, 25.0, 33.4]
TB = {'++': (1158.0, 1441.0), '+-': (1152.0, 1439.0)}
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)
def cfg(lam, tag, kf):
    LAM['v'] = lam; return DP.run(J, D, 'cv', [tag], kf, pr=lambda s: None)
W = {(l, kf): cfg(l, TAG, kf) for l in (0.0, LAM0, 1.0) for kf in (False, True)}
P = {kf: cfg(1.0, PTAG, kf) for kf in (False, True)}
pr(f"# dc_final (1100,100), record set {TAG}; lambda = {LAM0} +- {DLAM}; T_B (|Y|/|X|) ++ {TB['++']}, +- {TB['+-']}")
pr("| model | nodes | h | r (lambda-corrected) = r +- stat +- lambda | rms | L central (|Y|..|X|) | L at 2 sigma | P_J-cv twin: r +- stat, rms, L central, L 2 sigma | k=3 diff (+-): lambda-corr / P_J-cv |")
pr("|---|---|---|---|---|---|---|---|---|")
for model, kf, nus in (('F1 (A+C)', False, N4), ('F1 (A+C)', True, N6), ('F1B (A+B+C)', True, N6)):
    mk = 'F1B' if 'B' in model.split()[0] else 'F1'
    for h, sg in (('++', 1.0), ('+-', -1.0)):
        f = F.fit(W[(LAM0, kf)], h, nus, mk); f0 = F.fit(W[(0.0, kf)], h, nus, mk); f1 = F.fit(W[(1.0, kf)], h, nus, mk)
        sl = DLAM*abs(f1['r'] - f0['r']); st = np.hypot(f['er'], sl)
        a = abs(f['r']) + f['rms']
        fp = F.fit(P[kf], h, nus, mk); ap = abs(fp['r']) + fp['rms']
        k3 = []
        for o in (W[(LAM0, kf)], P[kf]):
            g = F.fit(o, h, nus, 'F1'); c3, c3j = c3fit(o, h, nus); R = c3j.shape[0]
            dj = c3j[:, 0] - sg*g['cj'][:, 0]
            k3.append(f"{c3[0] - sg*g['c'][0]:+.0f}+-{np.sqrt((R - 1)/R*np.sum((dj - dj.mean())**2)):.0f}")
        pr(f"| {model} | {len(nus)} | {h} | {f['r']:+.0f} +- {f['er']:.0f} +- {sl:.0f} (tot {st:.0f}) | {f['rms']:.0f} | "
           f"{TB[h][0]/a:.1f}-{TB[h][1]/a:.1f} | {TB[h][0]/(a + 2*st):.1f} | {fp['r']:+.0f} +- {fp['er']:.0f}, {fp['rms']:.0f}, "
           f"{TB[h][0]/ap:.1f}-{TB[h][1]/ap:.1f}, {TB[h][0]/(ap + 2*fp['er']):.1f} | " + (" / ".join(k3) if h == '+-' else '-') + " |")
open(os.path.join(dc_guard.CACHE, 'dc_final%s.txt' % TAG), 'w').write('\n'.join(lines) + '\n')
