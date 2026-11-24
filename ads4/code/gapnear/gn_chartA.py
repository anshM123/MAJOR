"""
gn_chartA (GAP_NEARGAP.md §3): chart-A test.  For every state (J, m/Lam) evaluated at Delta_g = 1000 (cache/gn_bind_base1000.json) and at
Delta_g = 2000 / 5000 (cache/gn_bind_<tag>.json, same protocol), compare at the SAME (J, m/Lam):
  delta = lam_K2 - lam_hybrid   (the beyond-leading-Regge correction, N_s units), the K2 step lam_K2 - lam_K1, and the budget terms,
and their ratios R = X(Dg)/X(1000), against the chart-A expectation (Lam_1000/Lam_Dg)^2 (0.2497 at 2000, 0.0399 at 5000).
Reads only gapnear/cache.  Usage: python3 gn_chartA.py tag2000[,tag...] > cache/chartA.txt
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
C = os.path.join(HERE, 'cache')
lam = lambda D: np.sqrt((D - 2.0)*(D - 1.0))


def load(tag):
    return json.load(open(os.path.join(C, 'gn_bind_%s.json' % tag)))


base = {(r['J'], round(r['mfl'], 5)): r for r in load('base1000')}
rows = []
for tag in sys.argv[1].split(','):
    rows += load(tag)
print("# chart A: same (J, m/Lam); delta = lam_K2 - lam_hyb; step = |lam_K2 - lam_K1|; ratios to Delta_g = 1000")
print(f"# expectation (Lam_1000/Lam)^2: 2000 -> {(lam(1000)/lam(2000))**2:.4f}, 5000 -> {(lam(1000)/lam(5000))**2:.4f}")
print(f"{'Dg':>5} {'J':>3} {'m/Lam':>8} | {'lam_K2':>9} {'delta':>10} {'R(delta)':>9} | {'step':>9} {'R(step)':>8} | "
      f"{'E_trunc':>9} {'R':>6} | {'E_low':>9} {'R':>6} | {'E_edge':>9} {'R':>6} | {'L':>7} {'L1000':>7} | {'L_lam':>7} {'L_lam1000':>9}")
for r in rows:
    k = (r['J'], round(r['mfl'], 5))
    b = base.get(k)
    d = r['lam_K2'] - r['lam_hybrid']; st = abs(r['lam_K2'] - r['lam_K1'])
    if b is None:
        print(f"{r['Dg']:>5} {r['J']:>3} {r['mfl']:8.5f} | {r['lam_K2']:+9.5f} {d:+10.2e} {'—':>9} | {st:9.2e} {'—':>8} | "
              f"{r['E_trunc']:9.2e} {'—':>6} | {r['E_low']:9.2e} {'—':>6} | {r['E_edge']:9.2e} {'—':>6} | {r['L']:7.1f} {'—':>7} | "
              f"{r['L_lam']:7.1f} {'—':>9}   (no 1000 state at this m/Lam: seam or new mass)")
        continue
    db = b['lam_K2'] - b['lam_hybrid']; sb = abs(b['lam_K2'] - b['lam_K1'])
    f = lambda a, c: (a/c if c else np.nan)
    print(f"{r['Dg']:>5} {r['J']:>3} {r['mfl']:8.5f} | {r['lam_K2']:+9.5f} {d:+10.2e} {f(d, db):9.3f} | {st:9.2e} {f(st, sb):8.3f} | "
          f"{r['E_trunc']:9.2e} {f(r['E_trunc'], b['E_trunc']):6.3f} | {r['E_low']:9.2e} {f(r['E_low'], b['E_low']):6.3f} | "
          f"{r['E_edge']:9.2e} {f(r['E_edge'], b['E_edge']):6.3f} | {r['L']:7.1f} {b['L']:7.1f} | {r['L_lam']:7.1f} {b['L_lam']:9.1f}")
