"""
e2b: independent check of the double-precision contour moments (e1, Delta_g = 100) by the Arb Laplace/Barnes method (skelexact3.py),
even k = 2, 4, 6, on all e1 states with m/Lam <= MFLMAX (default 1.5).  Prints max |d rho_k| (units of E_k) per state band and the
resulting change of the M02 action |dF|/Ns (even-k parts only).
Usage: python3 e2b_arbcheck.py Dg Df [MFLMAX]  -> cache/e2b_<Dg>_Df<Df>.txt
"""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
Dg = int(sys.argv[1]); DF = float(sys.argv[2]); MFLMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 1.5
os.environ['OMEGAU_DF'] = '%g' % DF
OMU = os.path.join(HERE, '..', 'omegau'); sys.path.insert(0, OMU)
import numpy as np
from common import load, gl01, peval, SIGMA
from adsaction import envk
from skelexact3 import Laplace3
F, PH = load('M02')
pn, pw = gl01(64)
Ek = lambda k, x: (2 - x*x)/(1 - x*x)**(1 + k/2)
db = pickle.load(open(os.path.join(HERE, 'cache', 'e1_%d_Df%g.pkl' % (Dg, DF)), 'rb'))
lines = []; worst = 0.0; worstF = 0.0
for key in sorted(db):
    if key[0] > MFLMAX + 1e-9:
        continue
    r = db[key]
    L3 = Laplace3(r['Delta'], key[1], DF, ks=(2, 4, 6), numax=float(r['nus'].max()))
    u = 1/key[0]**2; x = r['x']
    d = {}; dF = 0.0
    for k in (2, 4, 6):
        v = np.array([float(L3.value(k, nu).mid()) for nu in r['nus']])
        d[k] = (v - r['ex'][k])/Ek(k, x)
        dF += np.sum(u**(k - 2)*peval(PH[k], pn)*pw*envk(k, x)*d[k])
    md = max(np.abs(d[k]).max() for k in d)
    worst = max(worst, md); worstF = max(worstF, abs(dF))
    s = f"m/Lam={key[0]:.2f} J={key[1]:2d}: max_k |d rho| = {md:.1e};  |dF| (even k, M02 normal form) = {abs(dF):.1e}"
    lines.append(s)
lines.append(f"# worst |d rho| = {worst:.2e}, worst |dF| = {worstF:.2e}  ({len(lines)} states)")
print("\n".join(lines[-5:]))
open(os.path.join(HERE, 'cache', 'e2b_%d_Df%g.txt' % (Dg, DF)), 'w').write("\n".join(lines) + "\n")
