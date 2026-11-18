"""h13a_summary: route-A' hybrid tables and the two-term large-bhat envelope check (cache/h13a_hyb_*.pkl) -> cache/h13a_summary.txt.
Two-term form  A1 M1(beta)/Lam^2 + A2 M2(beta)/(2 Lam^3) = A1 M1/Lam^2 - A2/(2 bhat^3)  (M2/Lam^3 = -1/bhat^3 exactly),
A_n = 2 psi2^(n)(0) + sigma u psi3^(n)(0) (route A'; A1, A2 coincide with the B'-uplift values)."""
import os, sys, glob, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from math import factorial
from h13_common import CACHE, Lam_of, F
from scipy.interpolate import CubicSpline
BT, M1T = pickle.load(open(os.path.join(CACHE, 'h13_M12.pkl'), 'rb'))
cs = CubicSpline(np.log(BT), M1T); M1 = lambda b: cs(np.log(b))
d = lambda c, n: factorial(n)*float(c[n]) if n < len(c) else 0.0
An = lambda h, u, n: 2*d(F.s2, n) + (1 if h == '++' else -1)*u*d(F.s3, n)
out = []
def pr(s=''):
    print(s, flush=True); out.append(s)
pr(f"# route A': A1 = {d(F.s2,1)*2:.4f} + sigma u {d(F.s3,1):.4f};  A2 = {2*d(F.s2,2):.4f} + sigma u {d(F.s3,2):.3e};  A3 = {2*d(F.s2,3):.1f} + sigma u {d(F.s3,3):.1f}")
MB = ((1.5, 2.0), (2.0, 3.0), (3.0, 5.0), (5.0, 10.01), (10.01, 101))
BB = ((20, 40), (40, 70), (70, 120), (120, 200.01), (200, 400), (400, 1000.1))
for fn in sorted(glob.glob(os.path.join(CACHE, 'h13a_hyb_*_*_*.pkl'))):
    tag = os.path.basename(fn)[9:-4]; Dg = int(tag.split('_')[0]); L = Lam_of(Dg)
    recs = [r for v in pickle.load(open(fn, 'rb')).values() for r in v]
    pr(f"## Dg = {Dg}, bhat {tag.split('_')[1]}-{tag.split('_')[2]}: {len(recs)} records, {len(set(r['mfl'] for r in recs))} masses "
       f"[{min(r['mfl'] for r in recs):.3f}, {max(r['mfl'] for r in recs):.1f}]; non-positive {sum(r['Fh'] <= 0 for r in recs)};"
       f" min curv/Fh {min((r['curv2'] + r['curv3'])/r['Fh'] for r in recs):+.3e}")
    for h in ('++', '+-'):
        R = [r for r in recs if r['h'] == h]; w = min(R, key=lambda r: r['Fh']/r['Ns'])
        pr(f"   {h}: min F_hyb/N_s = {w['Fh']/w['Ns']:+.4f} at m/Lam={w['mfl']:.3f}, J={w['J']}, bhat={w['bh']:.1f}")
        for mb in MB:
            row = []
            for bb in BB:
                v = [r for r in R if mb[0] <= r['mfl'] < mb[1] and bb[0] <= r['bh'] < bb[1]]
                row.append(f"{min(r['Fh']/r['Ns'] for r in v):.3f}" if v else '  -  ')
            if any(c.strip() != '-' for c in row):
                pr(f"      m [{mb[0]:5.2f},{mb[1]:6.2f}): " + ' '.join(f"b[{bb[0]},{int(bb[1])}) {c}" for bb, c in zip(BB, row)))
    for b0 in (50, 100, 200):
        v = [r for r in recs if r['bh'] >= b0]
        if not v:
            continue
        rat = np.array([r['Fh']/(An(r['h'], 1/r['mfl']**2, 1)*M1(r['beta'])/L**2 - An(r['h'], 1/r['mfl']**2, 2)/(2*r['bh']**3)) for r in v])
        fl = np.array([r['Fflat']/(-An(r['h'], 1/r['mfl']**2, 2)/(2*r['bh']**3)) for r in v])
        pr(f"   bhat >= {b0}: F_hyb / two-term in [{rat.min():.5f}, {rat.max():.5f}];  F_flat/(-A2/(2 bhat^3)) in [{fl.min():.4f}, {fl.max():.4f}]  ({len(v)})")
open(os.path.join(CACHE, 'h13a_summary.txt'), 'w').write('\n'.join(out) + '\n')
