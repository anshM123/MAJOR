"""
dc_final2 (DIRCV.md §8): r and leverage at (1100, 100) with the new window nodes (nu = 11.0, 11.7, 12.3, 14.2, 15.0), record set TAG
(pooled), for: dir (lambda 0), cv (lambda 1), lambda-corrected (LAM0 +- DLAM), P_J-cv and P_J-dir (lambda-free).  Fits F1 (A+C) and
F1B (A+B+C) on (a) all pw window nodes nu in [11, 25], (b) the same + fw 16.7, 33.4, (c) (a) without nu = 13.3.
Also dr/d(eps) at nu = 13.3 (linear response, dir) to show how much the fit still hinges on that node.
Usage: DC_LOCAL=<tag>,<tag>P python3 dc_final2.py <tag> -> cache/dc_final2<tag>.txt
"""
import dc_guard, os, sys
import numpy as np
import z3_nodes as Z3
import dc_pipe as DP, dc_fit as F
LAM = {'v': 1.0}
_sc = Z3.scalar_check
def scalar_check(*a, **k):
    o = dict(_sc(*a, **k))
    o['cv'] = np.exp(LAM['v']*np.log(o['cv'])); o['cv_jk'] = np.exp(LAM['v']*np.log(o['cv_jk']))
    return o
Z3.scalar_check = scalar_check
J, D = 100, 1100.0
TAG = sys.argv[1]; PTAG = TAG + 'P'
LAM0, DLAM = float(os.environ.get('DC_LAM0', '0.6')), float(os.environ.get('DC_DLAM', '0.25'))
TB = {'++': (1158.0, 1441.0), '+-': (1152.0, 1439.0)}
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)
def cfg(lam, tag, kf, est='cv'):
    LAM['v'] = lam; return DP.run(J, D, est, [tag], kf, pr=lambda s: None)
W = {(l, kf): cfg(l, TAG, kf) for l in (0.0, LAM0, 1.0) for kf in (False, True)}
P = {kf: cfg(1.0, PTAG, kf) for kf in (False, True)}
Pd = {kf: cfg(1.0, PTAG, kf, est='dir') for kf in (False, True)}
nus_pw = [v for v in W[(0.0, False)]['nus'] if 10.5 <= v <= 25.01]
nus_fw = [v for v in W[(0.0, True)]['nus'] if 10.5 <= v <= 34.01]
sets = (('pw nodes', False, nus_pw), ('pw + fw', True, nus_fw), ('pw w/o 13.3', False, [v for v in nus_pw if abs(v - 13.3) > 1e-6]))
pr(f"# dc_final2 (1100,100) record set {TAG}; lambda = {LAM0} +- {DLAM};  pw nodes {np.round(nus_pw, 2).tolist()};  +fw {np.round(nus_fw, 2).tolist()}")
pr("| fit | nodes | h | dir | cv | lambda-corr r +- stat +- lam (rms) | L central / 2 sigma (lambda-corr) | P_J-cv r +- stat (rms) | L central / 2 sigma (P_J-cv) | P_J-dir r |")
pr("|---|---|---|---|---|---|---|---|---|---|")
for model in ('F1', 'F1B'):
    for name, kf, nus in sets:
        for h in ('++', '+-'):
            f0, f1, f = (F.fit(W[(l, kf)], h, nus, model) for l in (0.0, 1.0, LAM0))
            sl = DLAM*abs(f1['r'] - f0['r']); st = np.hypot(f['er'], sl); a = abs(f['r']) + f['rms']
            fp = F.fit(P[kf], h, nus, model); ap = abs(fp['r']) + fp['rms']; fd = F.fit(Pd[kf], h, nus, model)
            pr(f"| {model} | {name} ({len(nus)}) | {h} | {f0['r']:+.0f} +- {f0['er']:.0f} | {f1['r']:+.0f} +- {f1['er']:.0f} | "
               f"{f['r']:+.0f} +- {f['er']:.0f} +- {sl:.0f} ({f['rms']:.0f}) | {TB[h][0]/a:.1f}-{TB[h][1]/a:.1f} / {TB[h][0]/(a + 2*st):.1f} | "
               f"{fp['r']:+.0f} +- {fp['er']:.0f} ({fp['rms']:.0f}) | {TB[h][0]/ap:.1f}-{TB[h][1]/ap:.1f} / {TB[h][0]/(ap + 2*fp['er']):.1f} | "
               f"{fd['r']:+.0f} +- {fd['er']:.0f} |")
# chi2 of the fits (diagonal jackknife errors)
pr("")
for model in ('F1', 'F1B'):
    for name, kf, nus in sets[:2]:
        s = []
        for lab, o in (('lam', W[(LAM0, kf)]), ('PJcv', P[kf])):
            for h in ('++', '+-'):
                f = F.fit(o, h, nus, model); s.append(f"{lab} {h} chi2 {f['chi2']:.1f}/{f['dof']}")
        pr(f"  {model} {name}: " + "; ".join(s))
open(os.path.join(dc_guard.CACHE, 'dc_final2%s.txt' % TAG), 'w').write('\n'.join(lines) + '\n')
