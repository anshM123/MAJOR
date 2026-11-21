"""hc_sum: per-J summary of hc_chain outputs (w2/z4 modes: per-state records): min lambda/N_s, and lambda at the near-gap masses
m/Lam <= HC_MMAX (default 1.15), for one or more tags.  Usage: python3 hc_sum.py mode Dg tag1 [tag2 ...] [J=..]"""
import os, sys, pickle
import numpy as np
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache')
mode, Dg = sys.argv[1], int(sys.argv[2])
tags = [a for a in sys.argv[3:] if not a.startswith('J=')]
Js = [int(v) for a in sys.argv[3:] if a.startswith('J=') for v in a[2:].split(',')]
MM = float(os.environ.get('HC_MMAX', '1.15'))
for t in tags:
    o = pickle.load(open(os.path.join(C, 'hc_chain_%s_%s.pkl' % (mode, t)), 'rb'))[Dg]
    print(f"== {t}")
    for J in sorted(set(r['J'] for r in o)):
        if Js and J not in Js:
            continue
        R = [r for r in o if r['J'] == J]
        w = min(R, key=lambda r: r['lam'])
        ng = sorted((r['mfl'], r['lam'], r['diag']) for r in R if r['mfl'] <= MM + 1e-9)
        seen = {}
        for m, l, d in ng:
            k = round(m, 3)
            if k not in seen or l < seen[k][0]:
                seen[k] = (l, d)
        print(f"  J={J:2d}: min {w['lam']:+.4f} at {w['mfl']:.3f};  near gap: " +
              " ".join(f"{k:.3f}:{v[0]:+.4f}[" + "/".join(f"{x:+.3f}" for x in v[1]) + "]" for k, v in sorted(seen.items())))
