"""
lsw/lsw_classwise.py (LARGESPIN_W.md §4.2): class-wise tolerances K_tol[s] = H/c_s (one class at a time, sign-free, uniform in the
coupling ratio via H = f++ f+-/(f++ + f+-)), R1 and R2 forms, minimum over sampled masses per J; and the tolerated multiple of the
scalar-CHMRS coefficient vector (sign-free per class: sum_s |kappa_s^sc| c_s), plus the same with the nu^2 classes (B, C2) scaled by
beta^2 (the beyond-hybrid estimate of §3.4).
Usage: python3 lsw_classwise.py -> cache/classwise_1000.txt
"""
import json, os, sys
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'cache', 'costs_1000.json')))
SC = {'A ': 86.75/6 - 8.0, 'B ': 11/6, 'C1': 1/6, 'C2': 1/6, 'C3': 11.75/6}
out = []; pr = lambda s: (print(s, flush=True), out.append(s))
rows = {}
for r in R:
    f = {h: v['floor'] for h, v in r['fl'].items()}
    H = f['++']*f['+-']/(f['++'] + f['+-']) if len(f) == 2 else list(f.values())[0]
    c = r['cost']; b = r['beta']
    d = {}
    for k in ('A ', 'B ', 'C1', 'C2', 'C3', 'sh+1', 'sh-1', 'sh+2', 'sh-2'):
        d[k] = (H/abs(c[k][0]), H/abs(c[k][1]))
    for form in (0, 1):
        d['Lsc%d' % form] = H/sum(abs(SC[k])*abs(c[k][form]) for k in SC)
        d['Lb2%d' % form] = H/sum(abs(SC[k])*abs(c[k][form])*(b*b if k in ('B ', 'C2') else 1) for k in SC)
    rows.setdefault(r['J'], []).append(d)
pr('per J (min over sampled masses <= 1.3 Lam): K_tol R1/R2 for A, B, C1, C2, C3 ; index-shifted (R2) ; tolerated multiple of the scalar '
   'CHMRS vector R1/R2 ; same with nu^2 classes x beta^2 R1/R2')
for J in sorted(rows, key=lambda j: (j % 2, j)):
    m = {k: (min(x[k][0] for x in rows[J]), min(x[k][1] for x in rows[J])) if isinstance(rows[J][0][k], tuple) else min(x[k] for x in rows[J])
         for k in rows[J][0]}
    sh = min(m[k][1] for k in ('sh+1', 'sh-1', 'sh+2', 'sh-2'))
    pr('J=%3d | A %.2e/%.2e B %.2g/%.2g C1 %.2e/%.2e C2 %.2g/%.2g C3 %.2e/%.2e | sh %.2e | Lsc %.3g/%.3g | Lb2 %.3g/%.3g' %
       (J, *m['A '], *m['B '], *m['C1'], *m['C2'], *m['C3'], sh, m['Lsc0'], m['Lsc1'], m['Lb20'], m['Lb21']))
open(os.path.join(HERE, 'cache', 'classwise_1000.txt'), 'w').write('\n'.join(out) + '\n')
