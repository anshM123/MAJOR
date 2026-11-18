"""
certA_inv: the exact eta-inversion of psi_2 (spinfunc/h1_etainv_exact.py; read-only import of its pickle and of h2's helpers) on the states
of cert_actionA.py.  F_inv = F + dF,  dF = int_0^1 (psit2 - psi2) Mod2 dp  (h2_neargap_exact.py: flat + curvature part on h1's fine nu-grid,
remainder part on the 64 nodes with the e1 double-precision rho_2; tails p > 1 excluded [I]).  dF is [DV] (floating point; h1 checks
V1-V3), not certified.  Combined statement per state:  F_inv >= Fex_lo + dF - TOL |dF|  with TOL = 0.05 (allowance for h1/h2 numerics).
Usage: python3 certA_inv.py [Dg]  -> cache/certA_inv_<Dg>.txt, cache/certA_inv_<Dg>.pkl
"""
import os, sys, pickle, glob
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SPF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SPF)
OMU = os.path.join(HERE, '..', 'omegau'); sys.path.insert(0, OMU)
import numpy as np
from scipy.special import j0
import h2_neargap_exact as H2                       # Pcone, phit_at (read-only)
from common import gl01, Dh_one
from nf_common import NF, envk, SIGMA
TOL = 0.05
Dg = int(sys.argv[1]) if len(sys.argv) > 1 else 100
F = NF(); L = np.sqrt((Dg - 2.0)*(Dg - 1.0)); pn, pw = gl01(64)
d = pickle.load(open(os.path.join(SPF, 'cache', 'h1_%d.pkl' % Dg), 'rb'))
msk = d['nus'] <= L
nuf, wf = d['nus'][msk], d['wns'][msk]/L; pfine = nuf/L; dpf = (d['phit'] - d['phi'])[msk]
dpn = H2.phit_at(d, L*pn) - F.psi2(pn)
out, recs = [], []
def pr(s=''):
    print(s, flush=True); out.append(s)
pr(f"# certA_inv Dg={Dg}: exact eta-inversion shift dF on the certified states (TOL = {TOL})")
for Df in (3, 2):
    cfn = os.path.join(HERE, 'cache', 'certA_%d_Df%g.pkl' % (Dg, Df))
    if not os.path.exists(cfn):
        continue
    cert = pickle.load(open(cfn, 'rb'))
    e1 = {}
    for fn in glob.glob(os.path.join(HERE, 'cache', 'e1_%d_Df%g*.pkl' % (Dg, Df))):
        e1.update(pickle.load(open(fn, 'rb')))
    for (mfl, J), c in sorted(cert.items()):
        r0 = e1.get((round(mfl, 6), J)) or e1.get((mfl, J))
        if r0 is None:
            pr(f"  Df={Df} m/Lam={mfl} J={J}: no e1 record (rho) -> dF not computed"); continue
        u = 1/mfl**2; x = pn*np.sqrt(u); y = x*x; xf = pfine*np.sqrt(u); yf = xf*xf
        Pf = H2.Pcone(nuf, r0['eta'])
        for h, v in c['res'].items():
            fac = (1 - y)**2 if h == '++' else (1 - y)**-2
            Dh = Dh_one(h, J, yf); Dt = Dh if h == '++' else Dh/(1 - yf)**4
            T = (2 - yf)*Dt + 2*(Pf - j0((2*J + 1)*xf))
            dF = float(np.sum(wf*dpf*T) + np.sum(dpn*pw*envk(2, x)*fac*r0['rho'][2]))
            lo = v['Fex_lo'] + dF - TOL*abs(dF) if np.isfinite(v['Fex_lo']) else -np.inf
            recs.append(dict(Df=Df, mfl=mfl, J=J, h=h, Ns=v['Ns'], Fex_lo=v['Fex_lo'], F_lo=v['F_lo'], dF=dF, Finv_lo=lo, cert=v['cert']))
for Df in (3, 2):
    R = [r for r in recs if r['Df'] == Df]
    if not R:
        continue
    ok = [r for r in R if r['Finv_lo'] > 0]
    w = min((r for r in R if np.isfinite(r['Finv_lo'])), key=lambda r: r['Finv_lo']/r['Ns'])
    wd = max(R, key=lambda r: abs(r['dF'])/r['Ns'])
    pr(f"## Df={Df}: {len(R)} cases; Fex_lo - quad > 0 (certified, no inversion): {sum(r['cert'] for r in R)};"
       f"  with exact inversion (Fex_lo + dF - {TOL}|dF| > 0): {len(ok)}")
    pr(f"   min (Fex_lo + dF - TOL|dF|)/N_s = {w['Finv_lo']/w['Ns']:+.4f} [{w['h']} m/Lam={w['mfl']} J={w['J']}];"
       f"  max |dF|/N_s = {abs(wd['dF'])/wd['Ns']:.3e} [{wd['h']} m/Lam={wd['mfl']} J={wd['J']}]")
    for r in R:
        if not r['Finv_lo'] > 0:
            pr(f"   not established: {r['h']} m/Lam={r['mfl']} J={r['J']}: Fex_lo/Ns={r['Fex_lo']/r['Ns']:+.4f} dF/Ns={r['dF']/r['Ns']:+.4f}")
pickle.dump(recs, open(os.path.join(HERE, 'cache', 'certA_inv_%d.pkl' % Dg), 'wb'))
open(os.path.join(HERE, 'cache', 'certA_inv_%d.txt' % Dg), 'w').write('\n'.join(out) + '\n')
