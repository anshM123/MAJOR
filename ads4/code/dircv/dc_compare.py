"""dc_compare (DIRCV.md §3): r with jackknife errors for every estimator / record set at (1100, 100); paired same-CRN differences."""
import dc_guard, os, sys
import numpy as np
import dc_pipe as DP, dc_fit as F
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)

J, D = 100, 1100.0
cfg = {}
for est in ('cv', 'dir'):
    for tags, kf in ((['_w'], False), (['_w'], True), (['_wP'], False), (['_wP'], True), (['_none'], True)):
        cfg[(est, tags[0], kf)] = DP.get(J, D, est, tags, kf)
N4 = [13.3, 16.3, 20.8, 25.0]
pr("# (1100,100) C2, F1 (A+C) and F1B (A+B+C): r +- jackknife;  rms;  chi2/dof with diagonal jk errors")
for (est, tag, kf), o in cfg.items():
    nus = [v for v in o['nus'] if 12 <= v <= 34.0001]
    for model in ('F1', 'F1B'):
        for sel_name, sel in (('4 pw nodes', N4), ('all [12,34]', nus)):
            if tag == '_none' and sel_name == '4 pw nodes':
                continue
            if model == 'F1B' and len(sel) < 5:
                continue
            s = []
            for h in ('++', '+-'):
                f = F.fit(o, h, sel, model)
                s.append(f"{h}: r = {f['r']:+6.0f} +- {f['er']:4.0f} (c0 {f['c'][0]:+6.0f}), rms {f['rms']:4.0f}, chi2/dof {f['chi2']:.1f}/{f['dof']}")
            pr(f"  {est:3s} {tag:6s} keepfw={int(kf)} {model:3s} [{sel_name}: {np.round(sel, 1).tolist()}]  " + " | ".join(s))
pr("")
pr("# per-node y' = Gamma - kappa - (A-1)kappa  +- jk  (C2), 4 pw nodes")
for (est, tag, kf), o in cfg.items():
    if kf or tag == '_none':
        continue
    for h in ('++', '+-'):
        f = F.fit(o, h, N4)
        pr(f"  {est:3s} {tag:4s} {h}: " + ", ".join(f"{v:+6.0f}+-{e:4.0f}" for v, e in zip(f['y'], f['ey'])))
pr("")
pr("# paired same-CRN differences (jackknife of the difference), 4 pw nodes, F1:  Delta r;  Delta y per node")
for a, b in ((('cv', '_w', False), ('cv', '_wP', False)), (('dir', '_w', False), ('dir', '_wP', False)),
             (('dir', '_w', False), ('cv', '_wP', False)), (('cv', '_w', False), ('dir', '_wP', False)),
             (('dir', '_w', False), ('cv', '_w', False)), (('dir', '_wP', False), ('cv', '_wP', False))):
    for h in ('++', '+-'):
        p = F.paired(cfg[a], cfg[b], h, N4)
        chi = float(np.sum((p['dy']/p['edy'])**2))
        pr(f"  {a[0]}{a[1]} - {b[0]}{b[1]} {h}: Delta r = {p['dr']:+6.0f} +- {p['edr']:4.0f} ({abs(p['dr'])/p['edr']:.1f} sigma);  Delta y = "
           + ", ".join(f"{v:+6.0f}+-{e:4.0f}" for v, e in zip(p['dy'], p['edy'])) + f";  chi2(4) = {chi:.1f}")
open(os.path.join(dc_guard.CACHE, 'dc_compare.txt'), 'w').write('\n'.join(lines) + '\n')
