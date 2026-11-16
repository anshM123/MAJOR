"""Markdown tables for SPINNING_LOWK.md from cache/data_t1_l{1,2}.pkl: R_pres / R_flip with errors, k = 2, 3; m^2(R-1); and the
k = 4, 6 control (corner-completed vs e2e values).   Output: out_t3_tables.txt"""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np
from corner import CACHE
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
def fmt(v, e):
    return f"{v:.6f}({max(1, int(round(e*1e6)))})" if e < 5e-6 else (f"{v:.5f}({int(round(e*1e5))})" if e < 5e-5 else f"{v:.4f}({int(np.ceil(e*1e4))})")
for l in (1, 2):
    fn = os.path.join(CACHE, f'data_t1_l{l}.pkl')
    if not os.path.exists(fn): continue
    data = pickle.load(open(fn, 'rb'))
    for k in (2, 3):
        pr(f"### l = {l}, k = {k}: R_pres / R_flip (measured/predicted; error in units of the last digit)")
        for beta in (0.4, 1.0, 2.0):
            cs = sorted({key[3] for key in data if key[0] == beta and key[1] == k}, key=lambda c: c[0])
            if not cs: continue
            pr(f"| β ≈ {beta}, (Δ,J) → | " + " | ".join(f"({c[0]:.2f}, {c[1]})" for c in cs) + " |")
            pr("|---|" + "---|"*len(cs))
            pr("| m | " + " | ".join(f"{data[(beta, k, 0.0, c)]['m']:.1f}" for c in cs) + " |")
            for nu in (0.0, 1.0, 3.0):
                pr(f"| ν = {nu:.0f} | " + " | ".join(fmt(data[(beta, k, nu, c)]['Rp'], data[(beta, k, nu, c)]['dRp']) + " / "
                                                      + fmt(data[(beta, k, nu, c)]['Rf'], data[(beta, k, nu, c)]['dRf']) for c in cs) + " |")
            for nu in (0.0, 1.0, 3.0):
                pr(f"| m²(R−1), ν = {nu:.0f} | " + " | ".join(f"{data[(beta, k, nu, c)]['m']**2*(data[(beta, k, nu, c)]['Rp']-1):.2f} / "
                                                             f"{data[(beta, k, nu, c)]['m']**2*(data[(beta, k, nu, c)]['Rf']-1):.2f}" for c in cs) + " |")
            pr("")
    pr(f"### l = {l}: control k = 4, 6 and the old k = 2 values: max over cases of |R_new - R_old| (pres / flip), and of the propagated error of R_new")
    for k in (2, 3, 4, 6):
        dp = max(abs(v['Rp'] - v['Rp_old']) for key, v in data.items() if key[1] == k)
        df = max(abs(v['Rf'] - v['Rf_old']) for key, v in data.items() if key[1] == k)
        ep = max(v['dRp'] for key, v in data.items() if key[1] == k); ef = max(v['dRf'] for key, v in data.items() if key[1] == k)
        sh = max(max(v['lc']) for key, v in data.items() if key[1] == k)
        pr(f"  k={k}: max|new-old| {dp:.1e} / {df:.1e};  max err(new) {ep:.1e} / {ef:.1e};  max LC share of a mode {sh:.1e}")
    pr("")
open(os.path.join(HERE, 'out_t3_tables.txt'), 'w').write("\n".join(out) + "\n")
