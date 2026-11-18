"""
e4: p-quadrature check of the exact-skeleton action (M02): the 64-node Gauss rule of OMEGA_U versus 96 and 128 nodes, on
near-gap states.  Usage: python3 e4_quad.py Dg Df  -> cache/e4_<Dg>_Df<Df>.txt
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
Dg = int(sys.argv[1]); DF = float(sys.argv[2])
os.environ['OMEGAU_DF'] = '%g' % DF
OMU = os.path.join(HERE, '..', 'omegau'); sys.path.insert(0, OMU)
import numpy as np
from scipy.special import j0, eval_legendre
import remk
from common import load, state, Transforms, KS, gl01, flat_unstripped_parts, flat_parts
from adsaction import action
from skelexact import ExactState

Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
F, PH = load('M02'); TR = Transforms(PH); W0 = Transforms({k: np.array([1.0]) for k in KS})
pf, wf = gl01(400)
Ek = lambda k, x: (2 - x*x)/(1 - x*x)**(1 + k/2)
lines = []
for (mfl, J, h) in ((1.0, 0, '++'), (1.02, 2, '++'), (1.04, 4, '+-'), (1.08, 4, '+-'), (1.1, 0, '++'), (1.2, 2, '++'),
                    (1.2, 4, '+-'), (1.3, 2, '++'), (1.4, 2, '++'), (1.5, 0, '++')):
    st = state(Dg, Lam, mfl, J)
    if st is None:
        continue
    u = 1/mfl**2
    _, W, P0 = W0.regge(Lam, st['beta'])
    vals = []
    for NP in (64, 96, 128):
        pn, pw = gl01(NP)
        x = pn*np.sqrt(u)
        ex = ExactState(st['Delta'], J, DF).values(Lam*pn)
        Pd = remk.Pders_vec(Lam*pn, st['eta'], nmax=0)[0]
        sk = {}
        for k in KS:
            E = Ek(k, x)
            sk[k] = ((ex[k] - E*eval_legendre(J, 1 - 2*x*x) - 2*(Pd - j0((2*J + 1)*x)))/E, np.zeros(NP), np.ones(NP, bool))
        fw = remk.forward_remainders(J, st['m2'], st['eta'], W, P0)
        sk['fw'] = (fw[0], 0.0, fw[3], 0.0, True)
        Jsu, parts, FWu = flat_unstripped_parts(F, PH, h, u, J, pf, wf)
        Jss, B2, B3, FWs = flat_parts(F, h, u, J)
        i = int(np.searchsorted(Jsu, J)); i2 = int(np.searchsorted(Jss, J))
        r = action(F, PH, TR, h, st, J, mfl, Lam, pn, pw, sk, parts[i], FWu[i], (B2[i2], B3[i2], FWs[i2]))
        vals.append((r['Fh'] + r['S'])/r['Ns'])
    s = f"m/Lam={mfl:.2f} J={J} {h}: F_exact/Ns with NP = 64, 96, 128: " + ", ".join(f"{v:+.6f}" for v in vals)
    print(s, flush=True); lines.append(s)
open(os.path.join(HERE, 'cache', 'e4_%d_Df%g.txt' % (Dg, DF)), 'w').write("\n".join(lines) + "\n")
