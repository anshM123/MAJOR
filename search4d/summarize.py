"""Aggregate vaidya_scan outputs: python3 -I summarize.py runs/v_*.json"""
import json, sys
import numpy as np
names = ['MMI(3)', 'MMI(comp)', 'SA', 'n5#90', 'n5#72', 'n5#60a', 'n5#60b', 'n5#10']
worst = np.full(8, np.inf); kmin = np.inf; nvec = 0; mism = 0.0; oneside = 0; nocap = 0; nedge = 0; prop = 0
rows = []
for f in sys.argv[1:]:
    d = json.load(open(f))
    mm = [x if x is not None else np.inf for x in d['min_margin_per_orbit']]
    worst = np.minimum(worst, mm); kmin = min(kmin, d['K_min_slack']); nvec += d['n_vectors']
    mism = max(mism, d['two_end_mismatch'] if not d['one_sided'] else 0); oneside += len(d['one_sided'])
    nocap += sum(1 for x in d['ncap_branches'] if x == 0); nedge += len(d['ncap_branches']); prop += d.get('n_propagated', 0)
    rows.append((d['dlt'], d['t0'], d['kappa'], min(mm), d['K_min_slack'], len(d['one_sided']), d['n_vectors'],
                 max(d['ncap_branches']), round(d['seconds_total'])))
for r in sorted(rows):
    print('dlt=%-5s t0=%-4s kappa=%-5s minmargin=%+.2e Kmin=%+.2e onesided=%d nvec=%d maxcapbranches=%d %ds' % r)
print('runs', len(rows), 'vectors', nvec, 'edges', nedge, 'edges without cap', nocap, 'one-sided tubes', oneside,
      'propagated branches', prop)
print('max two-end area mismatch (pairs found from both ends):', mism)
print('min K-slack (SSA/WM instances):', kmin)
for n, w in zip(names, worst):
    print(f'  {n:10s} min margin {w:+.3e}')
