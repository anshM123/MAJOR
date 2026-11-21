"""Collect the reference Stage-0 functionals (from cache/lp_g4.jsonl) into cache/functionals_q.json (floating point, [NE]).
Convention: F = int psi2 B2 + int psi3 B3 - a (dB4|_0 - B5|_0) - c B4|_0, c = 1/2; bound g4 M^6/(8 pi G) <= A log(M/m_IR) + B - kappa |g3|^2 M^8."""
import json, os
from common_q import CACHE
L = [json.loads(l) for l in open(os.path.join(CACHE, 'lp_g4.jsonl'))]
want = {'Q0-full-robust': dict(d2=12, d3=12, r=3, delta=0.05, b3=True, matter=False, psi3zero=False),
        'Q0-min-robust (no B3)': dict(d2=12, d3=0, r=3, delta=0.05, b3=False, matter=False),
        'Q0-full-robust-matter': dict(d2=12, d3=12, r=3, delta=0.05, b3=True, matter=True),
        'Q0-full-optimum (delta=0, not robust)': dict(d2=12, d3=12, r=3, delta=0.0, b3=True, matter=False),
        'Q0-full-optimum-matter (delta=0)': dict(d2=12, d3=12, r=3, delta=0.0, b3=True, matter=True),
        'Q0-r5-psi3zero-robust (AdS-style weights)': dict(d2=8, d3=8, r=5, delta=0.05, b3=True, psi3zero=True, matter=False),
        'Q0-r5-min-robust (no B3)': dict(d2=8, d3=0, r=5, delta=0.05, b3=False, matter=False)}
out = {}
for name, sel in want.items():
    rows = [r for r in L if r.get('status') == 0 and not r.get('G3') and all(r.get(k) == v for k, v in sel.items())]
    r = rows[-1]
    out[name] = {k: r[k] for k in ['A', 'B', 'kappa', 'a', 'c', 'psi2', 'psi3', 'd2', 'd3', 'r', 'delta', 'matter', 'b3']}
    print('%-45s A=%.3f B=%.3f kappa=%.4f a=%.4f' % (name, r['A'], r['B'], r['kappa'], r['a']))
json.dump(out, open(os.path.join(CACHE, 'functionals_q.json'), 'w'), indent=1)
