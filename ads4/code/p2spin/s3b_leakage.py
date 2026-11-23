"""
p2spin/s3b_leakage.py (P2_SPINNING.md §3.3): size of the helicity leakage of the heavy symbols (modes outside the leading |h| = 2
classes, s2_main helicity frame), split into X,Y-independent constants and X,Y-dependent polynomial parts, per order and q.
Normalisation: relative to the structure's own leading amplitude. -> cache/s3b_leakage.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
import s3_spin as S3
OUT = []
pr = lambda *a: (print(' '.join(str(x) for x in a), flush=True), OUT.append(' '.join(str(x) for x in a)))
o = S3.load()
z = {s: 0 for s in S3.xs + S3.ys}
for qv, hb in o.items():
    pr(f'q = {qv} (beta = {float(2*sp.log(qv)):.4f})')
    for h, lead in (('same', (2, 2)), ('opp', (2, -2)), ('OO', (2, -2)), ('odd', (2, 2))):
        L0 = abs(complex(sp.N(hb[h][0][lead])))
        for od in (1, 2):
            cls = {'other |h|=2 class': [], '|k| < 2 modes': []}
            for md, c in hb[h][od].items():
                if md in S3.CLASSES['same'] + S3.CLASSES['opp']:
                    if md in hb[h][0]:
                        continue
                    key = 'other |h|=2 class'
                else:
                    key = '|k| < 2 modes'
                c0 = abs(complex(sp.N(c.subs(z))))/L0
                pol = sp.Poly(sp.expand(c - c.subs(z)), *S3.xs, *S3.ys)
                cx = max([abs(complex(sp.N(v))) for v in pol.coeffs()] or [0])/L0
                cls[key].append((c0, cx))
            pr(f'   {h:4s} order {od}: ' + '; '.join(f'{k}: max const {max([a for a, b in v] or [0]):.3g}, max X,Y-coef {max([b for a, b in v] or [0]):.3g}'
                                                  for k, v in cls.items()))
open(os.path.join(HERE, 'cache', 's3b_leakage.txt'), 'w').write('\n'.join(OUT) + '\n')
