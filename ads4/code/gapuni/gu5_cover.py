"""
gu5_cover (GAP_UNIFORM.md §3.4, Table C): sign-free chart-A growth allowance per spin.
  sup|Gamma|_lo, sup|Gamma|_hi : max over diagonal kinds/helicities, x in [0.05, 0.9], on the gapnear K2 rows at ~1000+J and ~2000+J
                                 (J = 7: no 1000 row, the 2007 row is used for both)
  G_cert, G_hyb               : minimum over the tabulated masses (twist edge, 1.02, 1.05, 1.1, 1.3) of the 1000 tolerance:
                                 1x1 (gu1_tol) for J = 0, 2 and odd J; Weyl 2x2 (gu1b_weyl) for even J >= 4
  F = G / max(sup_lo, sup_hi)  : allowed uniform growth factor of sup|Gamma| at EVERY Delta_g >= 1000 (Lemma T: tolerance ~ Lam^2)
Reads gapuni/cache (gu1_tol.json, gu1b_weyl.txt, gu2_growth.json). Usage: python3 gu5_cover.py -> cache/gu5_cover.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gu_guard  # noqa: F401
import glob
C = os.path.join(HERE, 'cache'); GN = os.path.join(HERE, '..', 'gapnear', 'cache')
H = {}
for fn in sorted(glob.glob(os.path.join(GN, 'gn_hmono_*.json'))):
    H.update(json.load(open(fn)))
tol = json.load(open(os.path.join(C, 'gu1_tol.json')))
W = {}
for l in open(os.path.join(C, 'gu1b_weyl.txt')).read().splitlines()[1:]:
    f = l.split()
    W.setdefault(int(f[0]), []).append((float(f[3]), float(f[7])))
out = [f"{'J':>3} | {'sup|G| ~1000':>12} {'sup|G| ~2000':>12} {'ratio':>6} | {'G_cert':>9} {'G_hyb':>9} | {'F_cert':>7} {'F_hyb':>7}"]
for J in sorted(H, key=int):
    Ji = int(J); sl = sh = 0.0
    for k, rows in H[J].items():
        if 'off' in k:
            continue
        for D, rr in rows.items():
            v = max(abs(r[1]) for r in rr if 0.05 - 1e-9 <= r[0] <= 0.9 + 1e-9)
            if 985 <= float(D) < 1100:
                sl = max(sl, v)
            elif 1900 < float(D) < 2100:
                sh = max(sh, v)
    if sl == 0.0:
        sl = sh
    if Ji % 2 == 0 and Ji >= 4:
        gc = min(a for a, b in W[Ji]); gh = min(b for a, b in W[Ji])
    else:
        rs = [t for t in tol if t['J'] == Ji]
        gc = min(t['1000']['Gcert'] for t in rs); gh = min(t['1000']['Ghyb'] for t in rs)
    base = max(sl, sh)
    out.append(f"{Ji:>3} | {sl:12.1f} {sh:12.1f} {sh/sl:6.3f} | {gc:9.3g} {gh:9.3g} | "
               f"{(gc/base if gc > 0 else float('nan')):7.2f} {gh/base:7.2f}")
print("\n".join(out))
open(os.path.join(C, 'gu5_cover.txt'), 'w').write("\n".join(out) + "\n")
