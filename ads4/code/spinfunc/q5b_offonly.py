"""q5b (S1_S4_FIXES.md §2.6): per-state min eigenvalue for selected J with the DIAGONAL graviton correction switched off (off-diagonal only)
-- used for J = 32 at Delta_g = 1000, whose K2 +- diagonal on the Delta = 1600 row breaks down (|Gamma_hh| ~ 1e6), while the off-diagonal is
regular.  Usage: Q3_JEVEN=28,32 Q3_EXTRA_M=... python3 q5b_offonly.py Dg -> stdout, cache/q5b_offonly_<Dg>.txt"""
import os, sys, io
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import q3_analysis as Q
Dg = int(sys.argv[1])
buf = io.StringIO()
def pr(s=""):
    print(s, flush=True); buf.write(s + "\n")
grids, out = Q.run([Dg], lambda s: None)
recs = [r for r in out[Dg]['recs'] if r['J'] in Q.JEVEN and r['J'] >= 20 and len(r['hs']) == 2]
for J in sorted(set(r['J'] for r in recs)):
    row = []
    for r in sorted([r for r in recs if r['J'] == J], key=lambda r: r['mfl']):
        r2 = dict(r, dF=np.array([[0.0, r['dF'][0, 1]], [r['dF'][1, 0], 0.0]]))
        row.append(f"{r['mfl']:.2f}:{Q.lam(r2, 1.0):+.3f}(full {Q.lam(r, 1.0):+.3g})")
    pr(f"J={J} [{'hyb' if any(r['DF'] == 'hyb' for r in recs if r['J'] == J) else 'e1'}] off-diagonal only: " + " ".join(row))
open(os.path.join(HERE, 'cache', 'q5b_offonly_%d.txt' % Dg), 'w').write(buf.getvalue())
