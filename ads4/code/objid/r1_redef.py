"""
r1: D-safe functional with the object identity made a definition at the scalar-skeleton level:
   Omega^{NF,C} = int dp sum_k phi_k(p) eps_k Lam^{2k} C_{k,Lam p}   (CHMRS dispersive C_{k,nu}; spinning extension for gravitons, see
   FUNCTIONAL_REDEFINITION.md), no forward terms.
On the scalar skeleton the even-k terms ARE exactskel's exact C_{k,nu} (no model).  Odd k has no scalar functional; two readings bracket it:
   B' : odd density sigma (1-x^2)^{-(k+1)/2} with the continued scalar remainder (= exactskel e7, F_ex);
   E  : odd moments = eps_k x the exact k-continued scalar C (shift of e3's dodd_E, here with NF02 weights; 64 nodes).
Thresholds n(delta) and bounds (e9 formula, S_ff residue factor) at Delta_g = 100, 300, 1000 on all e1 states (both Df, high-J samples).
Usage: python3 r1_redef.py  -> cache/r1_redef.txt
"""
import os, sys, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ('omegau', 'exactskel'):
    sys.path.insert(0, os.path.join(HERE, '..', p))
import numpy as np
from scipy.integrate import quad
from common import gl01, Dh_one, EPS, SIGMA, KS
from nf_common import NF
import evenkernel as EK
pn, pw = gl01(64)
F = NF()
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

def dodd_E(h, J, mfl, rho):
    s = SIGMA[h]; u = 1/mfl**2; y = u*pn*pn; D = Dh_one(h, J, y); tot = 0.0
    for k in KS:
        if k % 2 == 0: continue
        Ek = (2 - y)/(1 - y)**(1 + k/2)
        tot += u**(k - 2)*np.sum(pw*F.phi(k, pn)*s*(EPS[k]*Ek - 1/(1 - y)**((k + 1)/2))*(D + rho[k]))
    return tot

def bound(Dg, n):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); L = Lam/n
    LA = quad(lambda t: float(EK.residue_factor_C2(Lam*t, 'Sff'))*float(F.psi2(t))/(t*t + 2.25/L**2), 0, 1, limit=400, points=[1.5/L])[0]
    return 240*np.sqrt(n**8*16*LA)/Lam**4

pr("# r1: Omega^{NF,C} (CHMRS C_{k,nu} by definition) on the exact scalar skeleton; odd k bracketed by readings B' and E")
for Dg in (100, 300, 1000):
    e7 = pickle.load(open(os.path.join(HERE, '..', 'exactskel', 'cache', 'e7_%d_recs.pkl' % Dg), 'rb'))
    e1 = {}
    for fn in glob.glob(os.path.join(HERE, '..', 'exactskel', 'cache', 'e1_%d_Df*.pkl' % Dg)):
        DF = float(os.path.basename(fn).split('_Df')[1].replace('.pkl', '').replace('_hiJ', ''))
        for key, r0 in pickle.load(open(fn, 'rb')).items():
            e1[(DF,) + tuple(key)] = r0
    recs = []
    for r in e7:
        r0 = e1.get((float(r['Df']), r['mfl'], r['J']))
        if r0 is None: continue
        dE = dodd_E(r['h'], r['J'], r['mfl'], r0['rho'])
        recs.append(dict(r, FB=r['Fex'], FE=r['Fex'] + dE))
    grid = sorted(set(r['mfl'] for r in recs))
    pr(f"## Delta_g = {Dg}: {len(recs)} (state, helicity, skeleton) records, m/Lam in [{grid[0]}, {grid[-1]}]")
    for rd in ('FB', 'FE'):
        for delta in (0.0, 0.03):
            n = next((n0 for n0 in grid if all(r[rd] - delta*r['Ns'] > 0 for r in recs if r['mfl'] >= n0 - 1e-9)), None)
            w = min(recs, key=lambda r: r[rd]/r['Ns'])
            pr(f"   reading {rd[1]} delta={delta}: n = {n}; min F/N_s = {w[rd]/w['Ns']:+.4f} [{w['h']} m/Lam={w['mfl']} J={w['J']} Df={w['Df']}]"
               + (f"; sqrt(t4^2+aT^2) <= {bound(Dg, n):.3g}" if n else ""))
open(os.path.join(HERE, 'cache', 'r1_redef.txt'), 'w').write("\n".join(out) + "\n")
