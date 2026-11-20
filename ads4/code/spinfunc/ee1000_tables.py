"""
ee1000_tables (EE_DG1000_LARGEJ.md): per-J tables from ee1000_chain runs (cache/ee1000_recs_<prefix>_<variant>.pkl): min lambda/N_s and
max |dF_hh|/N_h per variant, the ensemble minimum, and the per-state map (main variant).
Usage: python3 ee1000_tables.py prefix variant1 variant2 ...   -> cache/ee1000_tables_<prefix>.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
C = os.path.join(HERE, 'cache')
P, VS = sys.argv[1], sys.argv[2:]
R = {}
for v in VS:
    fn = os.path.join(C, 'ee1000_recs_%s_%s.pkl' % (P, v))
    if os.path.exists(fn):
        R[v] = pickle.load(open(fn, 'rb'))
lines = []
def pr(s=''):
    print(s); lines.append(s)
vs = [v for v in VS if v in R]
Js = sorted(set(r['J'] for v in vs for r in R[v]['recs']))
pr(f"# per J: min lam/N_s over states [m/Lam of the minimum] and max |dF_hh|/N_h, per variant; variants: {vs}")
pr("| J | uncorr. | " + " | ".join(vs) + " | ensemble min | max |dF|/N (all variants) |")
pr("|---|---|" + "---|"*len(vs) + "---|---|")
for J in Js:
    cells, mins, dmax = [], [], 0.0
    lam0 = None
    for v in vs:
        rj = [r for r in R[v]['recs'] if r['J'] == J]
        if not rj:
            cells.append('—'); continue
        w = min(rj, key=lambda r: r['lam'])
        d = max(abs(x) for r in rj for x in r['diag'])
        lam0 = min(r['lam0'] for r in rj)
        cells.append(f"`{w['lam']:+.3f}` [{w['mfl']:.3f}] ({d:.3f})")
        mins.append(w['lam']); dmax = max(dmax, d)
    pr(f"| {J} | `{lam0:+.3f}` | " + " | ".join(cells) + f" | **`{min(mins):+.3f}`** | {dmax:.3f} |")
if 'main' in R:
    pr("\n# per-state map, main variant: lam/N_s (diag corrections dF_hh/N_h)")
    for J in Js:
        rj = [r for r in R['main']['recs'] if r['J'] == J]
        pr(f"  J={J}: " + "  ".join(f"{r['mfl']:.3f}:{r['lam']:+.3f}(" + "/".join(f"{d:+.3f}" for d in r['diag']) + ")" for r in rj))
    pr("\n# rows used (main): " + "; ".join(f"J={J}: " + ", ".join(f"{D:.0f} ({i['n']} nodes, x {i['xs'][0]:.3f}-{i['xs'][-1]:.3f})" for D, i in R['main']['info'][J].items())
                                       for J in sorted(R['main']['info'])))
open(os.path.join(C, 'ee1000_tables_%s.txt' % P), 'w').write("\n".join(lines) + "\n")
