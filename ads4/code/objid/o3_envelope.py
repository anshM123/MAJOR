"""
o3: stress test of H14 beyond the three candidate objects: bulk-point ratio families R(k,x) = (1-r*^2)^alpha (all k, or k >= 3 only),
alpha in [-6, 12], on the NF02 (D-safe) and M02 (D-strong) exact-skeleton near-gap actions (model of o2; Delta_g = 100, both Df, all e1 states).
For each alpha: n(delta = 0) and min F/N_s.  alpha = 0 is the CHMRS C baseline; Pi ~ k, Psi ~ k+2 (times 1/(1+r^2)).
Usage: python3 o3_envelope.py -> cache/o3_envelope.txt
"""
import os, sys, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import o2_ads as O2
from o1_flat import rstar
from common import KS
ALPHAS = (-6, -4, -2, 2, 4, 8, 12)
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)
e7 = pickle.load(open(os.path.join(HERE, '..', 'exactskel', 'cache', 'e7_100_recs.pkl'), 'rb'))
e3 = pickle.load(open(os.path.join(HERE, '..', 'exactskel', 'cache', 'e3_100_recs.pkl'), 'rb'))
recs = []
for fn in sorted(glob.glob(os.path.join(HERE, '..', 'exactskel', 'cache', 'e1_100_Df*.pkl'))):
    DF = float(os.path.basename(fn).split('_Df')[1].replace('.pkl', '').replace('_hiJ', ''))
    for (mfl, J), r0 in pickle.load(open(fn, 'rb')).items():
        for h in ('++', '+-'):
            if (h == '++' and J % 2) or (h == '+-' and J < 4):
                continue
            for fun in ('NF02', 'M02'):
                ref = ([r for r in e7 if r['h'] == h and abs(r['mfl'] - mfl) < 1e-9 and r['J'] == J and float(r['Df']) == DF] if fun == 'NF02'
                       else [r for r in e3.get(DF, {}).get(('M02', h), []) if abs(r['mfl'] - mfl) < 1e-9 and r['J'] == J])
                if not ref:
                    continue
                x, I = O2.integrand(fun, h, r0); r2 = rstar(x)**2
                rec = dict(fun=fun, h=h, mfl=mfl, J=J, Ns=ref[0]['Ns'], Fex=ref[0]['Fex'])
                for a in ALPHAS:
                    rec[('all', a)] = sum((((1 - r2)**a - 1)*I[k]).sum() for k in KS)
                    rec[('k>=3', a)] = sum((((1 - r2)**a - 1)*I[k]).sum() for k in KS if k >= 3)
                recs.append(rec)
grid = sorted(set(r['mfl'] for r in recs))
pr("# o3: bulk-point-ratio family R = (1-r*^2)^alpha; n(delta=0), min F/N_s [worst state]")
for fun in ('NF02', 'M02'):
    for which in ('all', 'k>=3'):
        for a in ALPHAS:
            v = [r for r in recs if r['fun'] == fun]
            F = lambda r: r['Fex'] + r[(which, a)]
            n = next((n0 for n0 in grid if all(F(r) > 0 for r in v if r['mfl'] >= n0 - 1e-9)), None)
            w = min(v, key=lambda r: F(r)/r['Ns'])
            pr(f"  {fun} {which:5s} alpha={a:+d}: n = {n}; min F/N_s = {F(w)/w['Ns']:+.4f} [{w['h']} m/Lam={w['mfl']} J={w['J']}]")
open(os.path.join(HERE, 'cache', 'o3_envelope.txt'), 'w').write("\n".join(out) + "\n")
