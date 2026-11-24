"""
gn_hsum (GAP_NEARGAP.md §3.2): compact H_mono-A summary from cache/gn_hmono_*.json: per J, over all entries and x in [0.05, 0.9]:
max|Gamma| on the 2000 (5000) row, max|Gamma_2000 - Gamma_1000|, max|Gamma_2000 - Gamma_1600| (c2), and the K2 - K1 step on both rows.
Usage: python3 gn_hsum.py > cache/hsum.txt
"""
import os, glob, json
import numpy as np
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache')
H = {}
for fn in sorted(glob.glob(os.path.join(C, 'gn_hmono_*.json')), key=os.path.getmtime):   # newest last
    H.update(json.load(open(fn)))
print(f"{'J':>3} {'rows':>24} | {'max|G_hi|':>9} {'max|G_hi-G_1000|':>16} {'(x<=0.8)':>9} {'max|G_hi-G_1600|':>16} | {'step_1000':>9} {'step_hi':>8}")
for J in sorted(H, key=int):
    keys = H[J]
    Ds = sorted({D for k in keys for D in keys[k]}, key=float)
    lo = [D for D in Ds if 985 <= float(D) < 1100]
    mid = [D for D in Ds if float(D) == 1600]
    for hi in [D for D in Ds if float(D) > 1900]:
        a = b = c = s0 = s1 = b8 = 0.0
        for k, rows in keys.items():
            if 'off' in k:
                continue
            vh = {round(x, 4): (v, s) for x, v, s in rows.get(hi, [])}
            for x, (v, s) in vh.items():
                a = max(a, abs(v)); s1 = max(s1, s)
            for ref, slot in ((lo, 'lo'), (mid, 'mid')):
                if not ref:
                    continue
                for x, v, s in rows.get(ref[0], []):
                    if round(x, 4) in vh:
                        d = abs(vh[round(x, 4)][0] - v)
                        if slot == 'lo':
                            b = max(b, d); s0 = max(s0, s)
                            if x <= 0.8 + 1e-9:
                                b8 = max(b8, d)
                        else:
                            c = max(c, d)
        print(f"{J:>3} {('/'.join(lo + mid + [hi])):>24} | {a:9.1f} {b if lo else float('nan'):16.1f} {b8 if lo else float('nan'):9.1f} {c if mid else float('nan'):16.1f} | "
              f"{s0 if lo else float('nan'):9.0f} {s1:8.0f}")
