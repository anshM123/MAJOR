"""
e7: Omega^{NF,SC} (NF02 uplift, no forward terms; OMEGA_NF.md) with the exact scalar skeleton on all e1 states (e1_<Dg>_Df*.pkl, including
the high-J file *_hiJ.pkl).  Same model as omegau/nf_common.action_nf with rho_k exact on every node (no control cut).
Prints min F_exact/Ns per (m/Lam, bhat) band, the old series kappa on the same states, and thresholds n_exact(delta) restricted to
bhat <= 10, 20, 40 (only the scanned states: J <= 10 everywhere, J up to 28 on the high-J sample).
Usage: python3 e7_nf.py Dg  -> cache/e7_<Dg>.txt
"""
import os, sys, pickle, glob
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OMU = os.path.join(HERE, '..', 'omegau'); sys.path.insert(0, OMU)
os.environ.setdefault('OMEGAU_DF', '3')
import numpy as np
from nf_common import NF, TransformsNF, flat_parts_nf, action_nf, KS, SIGMA
from common import gl01

Dg = int(sys.argv[1])
pn, pw = gl01(64)
Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
F = NF(); TR = TransformsNF(F)
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

recs = []
for fn in sorted(glob.glob(os.path.join(HERE, 'cache', 'e1_%d_Df*.pkl' % Dg))):
    DF = os.path.basename(fn).split('_Df')[1].replace('.pkl', '').replace('_hiJ', '')
    db = pickle.load(open(fn, 'rb'))
    for (mfl, J), r0 in db.items():
        u = 1/mfl**2
        st = dict(m=r0['m'], m2=r0['m2'], eta=r0['eta'], beta=r0['beta'], Delta=r0['Delta'])
        for h in SIGMA:
            if (h == '++' and J % 2) or (h == '+-' and J < 4):
                continue
            pf, wf = gl01(max(400, 2*J + 100))
            Js, B2, B3, parts = flat_parts_nf(F, h, u, J, pf, wf)
            i = int(np.searchsorted(Js, J))
            ske = {k: (r0['rho'][k], np.zeros(64), np.ones(64, bool)) for k in KS}
            re = action_nf(F, TR, h, st, J, mfl, Lam, pn, pw, ske, parts[i], B2[i], B3[i])
            rs = action_nf(F, TR, h, st, J, mfl, Lam, pn, pw, r0['ser'], parts[i], B2[i], B3[i])
            recs.append(dict(Df=DF, h=h, mfl=mfl, J=J, bh=(2*J + 1)/mfl, Ns=re['Ns'], Fh=re['Fh'], Fex=re['Fh'] + re['S'],
                             kap=(rs['Fh'] - rs['E'])/rs['U'] if rs['U'] > 0 else np.inf, E=rs['E'], U=rs['U']))
pickle.dump(recs, open(os.path.join(HERE, 'cache', 'e7_%d_recs.pkl' % Dg), 'wb'))
pr(f"# e7: NF02 / Omega^NF,SC with the exact skeleton, Delta_g = {Dg}; {len(recs)} (state, helicity, skeleton) records")
pr("  band (m/Lam) x (bhat): n, min F_exact/Ns [state], min F_hyb/Ns, min series kappa")
for mb in ((1.0, 1.1), (1.1, 1.3), (1.3, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.01)):
    for bb in ((0, 5), (5, 10), (10, 15), (15, 20), (20, 30), (30, 60)):
        v = [r for r in recs if mb[0] <= r['mfl'] < mb[1] and bb[0] < r['bh'] <= bb[1]]
        if not v:
            continue
        w = min(v, key=lambda r: r['Fex']/r['Ns'])
        pr(f"   m[{mb[0]:.1f},{mb[1]:.1f}) bhat({bb[0]},{bb[1]}]: {len(v):4d}  {w['Fex']/w['Ns']:+.4f} [{w['h']} m/Lam={w['mfl']:.2f} J={w['J']} Df={w['Df']}]"
           f"  {min(r['Fh']/r['Ns'] for r in v):+.4f}  {min(r['kap'] for r in v):.3g}")
grid = sorted(set(r['mfl'] for r in recs))
for bmax in (5, 10, 20, 40, 1e9):
    for delta in (0.0, 0.03):
        v = [r for r in recs if r['bh'] <= bmax]
        n = next((n0 for n0 in grid if all(r['Fex'] - delta*r['Ns'] > 0 for r in v if r['mfl'] >= n0 - 1e-9)), None)
        bad = [r for r in v if r['Fex'] - delta*r['Ns'] <= 0]
        pr(f"  bhat <= {bmax:g}, delta = {delta}: n_exact = {n}  ({len(bad)} non-positive records" +
           (f"; highest m/Lam {max(r['mfl'] for r in bad):.2f}, bhat range {min(r['bh'] for r in bad):.1f}-{max(r['bh'] for r in bad):.1f})" if bad else ")"))
open(os.path.join(HERE, 'cache', 'e7_%d.txt' % Dg), 'w').write("\n".join(out) + "\n")
