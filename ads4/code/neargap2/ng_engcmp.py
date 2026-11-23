"""
ng_engcmp (NEARGAP2.md §1.3): the binding state's small-nu window.  Parity-odd analytic-radial engine (EE2, ng_engine5.py) against the
own-mass K2 row (and the direct points where they overlap), in the convention of po_k2cmp (PARITY_ODD_CLOSE §2.1):
   m^2 [D(nu) - D(nu0)],  D = log(state/reference),  nu0 = 0.05 m  (normalisation-free).
Engine truncation: Kt = 1, 2, 3 (differences).  Output: cache/ng_engcmp.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, SF); sys.path.insert(0, HERE)
import ng_guard  # noqa: F401
os.chdir(SF)
import numpy as np
import bg_eeq as EQ
import p4_odd_analysis as P4
NGC = os.path.join(HERE, 'cache')
J, D = 5, 1005.0
m = np.sqrt((D - J - 2)*(D + J - 1))
st = pickle.load(open(os.path.join(NGC, 'ng_po_eeq2_D1005_J5_K3.pkl'), 'rb'))
rf = pickle.load(open(os.path.join(NGC, 'ng_po_eeq2_D20000_J5_K3.pkl'), 'rb'))
nus = st['nus']
Deng = {Kt: np.log(EQ.value(st, 'EE2', 0, 0, Kt)[0]/EQ.value(rf, 'EE2', 0, 0, Kt)[0]) for Kt in (1, 2, 3)}
R = {(k[0], round(k[1], 9)): v for k, v in pickle.load(open(os.path.join(NGC, 'w1_p3_EE_J5_ng.pkl'), 'rb')).items() if k[2] == 1.0}


def _get(DD, nu):
    c = [k for k in R if k[0] == DD and abs(k[1] - nu) < 2e-6*max(nu, 1)]
    return R[c[0]] if c else None


def dk2(nu, K):
    s = _get(D, nu); r = _get(20000.0, nu)
    if s is None or r is None:
        return None
    dl = s['logC'] - s['logK'] - r['logC'] + r['logK']
    o = ('M0', 'M2', 'M4')
    return dl + np.log(P4.pred([s[k][-1, -1] for k in o], [r[k][-1, -1] for k in o], K))


i0 = int(np.argmin(np.abs(nus - 0.05*m)))
nu0 = nus[i0]
L = [f"# ng_engcmp J=5 Delta=1005 (m={m:.3f}): m^2[D(nu) - D(nu0)], nu0 = {nu0:.3f} (x = 0.05); engine Kt=3 [Kt=3 - Kt=2] | K2 (K=2) [K2 - K1]"]
k2ref = {K: dk2(nu0, K) for K in (1, 2)}
for i, nu in enumerate(nus):
    e3 = m*m*(Deng[3][i] - Deng[3][i0]); e2 = m*m*(Deng[2][i] - Deng[2][i0])
    s = f"  nu={nu:8.3f} x={nu/m:.4f}: engine {e3.real:+10.1f} [{(e3 - e2).real:+8.1f}]"
    v2 = dk2(nu, 2)
    if v2 is not None and k2ref[2] is not None:
        a2 = m*m*(v2 - k2ref[2]); a1 = m*m*(dk2(nu, 1) - k2ref[1])
        s += f" | K2 {a2.real:+10.1f} [{(a2 - a1).real:+8.1f}] | engine - K2 {(e3 - a2).real:+8.1f}"
    L.append(s)
print("\n".join(L))
open(os.path.join(NGC, 'ng_engcmp.txt'), 'w').write("\n".join(L) + "\n")
