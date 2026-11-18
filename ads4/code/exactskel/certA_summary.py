"""certA_summary: tables for CERTIFIED_NEARGAP_A.md from cache/certA_100_Df*.pkl and cache/certA_inv_100.pkl -> cache/certA_100_summary.txt"""
import os, sys, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
inv = {(r['Df'], r['mfl'], r['J'], r['h']): r for r in pickle.load(open(os.path.join(HERE, 'cache', 'certA_inv_100.pkl'), 'rb'))}
out = []
def pr(s=''):
    print(s, flush=True); out.append(s)
BANDS = ((1.01, 1.05), (1.05, 1.1), (1.1, 1.2), (1.2, 1.5), (1.5, 2.01), (2.4, 3.01))
for Df in (3, 2):
    d = pickle.load(open(os.path.join(HERE, 'cache', 'certA_100_Df%d.pkl' % Df), 'rb'))
    R = [(k, h, v) for k, r in d.items() for h, v in r['res'].items()]
    lowJ = [a for a in R if a[0][1] <= 10]; hiJ = [a for a in R if a[0][1] > 10]
    pr(f"## Df={Df}: {len(d)} states, {len(R)} (state, channel) cases (J <= 10: {len(lowJ)}, J = 11-16: {len(hiJ)})")
    pr(f"   F64 > 0 certified: {sum(v['cert64'] for _, _, v in R)};  exact integral F > 0 certified (F64_lo - quad.err > 0): {sum(v['cert'] for _, _, v in R)}")
    fin = [a for a in R if np.isfinite(a[2]['Fex_lo'])]
    w = min(fin, key=lambda a: a[2]['Fex_lo']/a[2]['Ns'])
    wm = min(R, key=lambda a: a[2]['F_mid']/a[2]['Ns'])
    pr(f"   min certified lower bound Fex_lo/N_s = {w[2]['Fex_lo']/w[2]['Ns']:+.4f} [{w[1]} m/Lam={w[0][0]} J={w[0][1]}];"
       f"  min midpoint F/N_s = {wm[2]['F_mid']/wm[2]['Ns']:+.4f} [{wm[1]} m/Lam={wm[0][0]} J={wm[0][1]}]")
    for b in BANDS:
        v = [a for a in fin if b[0] <= a[0][0] < b[1] and a[0][1] <= 10]
        if not v:
            continue
        w = min(v, key=lambda a: a[2]['Fex_lo']/a[2]['Ns'])
        e = max(a[2]['err']/a[2]['Ns'] for a in v)
        wi = min(v, key=lambda a: inv[(Df, a[0][0], a[0][1], a[1])]['Finv_lo']/a[2]['Ns'])
        ri = inv[(Df, wi[0][0], wi[0][1], wi[1])]
        pr(f"   m/Lam [{b[0]:.2f},{b[1]:.2f}): {len(v):3d} cases; min Fex_lo/N_s {w[2]['Fex_lo']/w[2]['Ns']:+.4f} [{w[1]} J={w[0][1]} m={w[0][0]}];"
           f" max quad.err/N_s {e:.1e};  with exact inversion (Fex_lo + dF - 0.05|dF|)/N_s >= {ri['Finv_lo']/ri['Ns']:+.4f} [{wi[1]} J={wi[0][1]} m={wi[0][0]}]")
    if hiJ:
        w = min(hiJ, key=lambda a: a[2]['Fex_lo']/a[2]['Ns'])
        pr(f"   J = 11-16: all {sum(v['cert'] for _, _, v in hiJ)}/{len(hiJ)} certified; min Fex_lo/N_s {w[2]['Fex_lo']/w[2]['Ns']:+.4f} [{w[1]} J={w[0][1]} m={w[0][0]}]")
    rad = max(((v['F_hi'] - v['F_lo'])/2/v['Ns'], k, h) for k, h, v in R)
    pr(f"   max ball radius / N_s: {rad[0]:.3f} [{rad[2]} m/Lam={rad[1][0]} J={rad[1][1]}]")
open(os.path.join(HERE, 'cache', 'certA_100_summary.txt'), 'w').write('\n'.join(out) + '\n')
