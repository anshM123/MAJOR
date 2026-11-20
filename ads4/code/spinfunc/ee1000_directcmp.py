"""
ee1000_directcmp (EE_DG1000_LARGEJ.md §3): direct (RQMC, ee1000_direct) STATE values vs K2 Laplace values of Gamma_{EE,h} at the same
(Delta, nu), both with the K2 heavy reference (the n1 contour masks most samples at the reference for J >~ 36: w2 W2_DSTATE convention,
EE_NEARGAP §5), the EE-own helicity basis, the Jacobi fix and the ee1000 loaders.  The direct state record replaces M0 by the all-orders
matrix (TPS normalisation) and sets M2 = M4 = 0 (w2._dstate).  Gamma of the direct variant is quoted with the K2 row's a_h.
Usage: python3 ee1000_directcmp.py "J:Delta,J:Delta"   -> cache/ee1000_directcmp.txt, .pkl
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import ee1000_chain as CH
import ee1000_direct as ED
import n1_analysis as N1
Q, W, C = CH.Q, CH.W, CH.C
for fn in sorted(glob.glob(os.path.join(C, '*shadow_k2*.pkl'))):
    N1.SHK2.update(pickle.load(open(fn, 'rb')))
USED = []


def _dstate(R, J, kind):
    if not W.DSTATE:
        return R
    db = ED.dload(J, kind)
    R = dict(R)
    I2 = dict(Ti=np.eye(2 if J >= 4 else 1))
    for (D, nu), r in db.items():
        if D >= 1e4:
            continue
        k = (float(D), round(float(nu), 9))
        if k not in R:
            continue
        H, l = N1.Hrep(r, I2, D, J)
        st = R[k]
        M = np.exp(l - (st['logC'] - st['logK']))*H.mean(0)
        Hs = np.exp(l - (st['logC'] - st['logK']))*H
        R[k] = dict(st, M0=M, M2=0*M, M4=0*M, M0reps=Hs)
        USED.append((kind, J, D, round(nu, 2)))
    return R


W._dstate = _dstate
_q2o = W._q2_even_orig                # DSTATE branch of w2._load reads C3H through this name: include the ee1000 caches
W._q2_even_orig = lambda kind, J: {**_q2o(kind, J), **W._pk(['ee1000_pts_%s_J%d_k2*.pkl' % (kind, J)])}
lines = []
def pr(s=""):
    print(s, flush=True); lines.append(s)
pr("# ee1000_directcmp: direct (RQMC) state / K2 reference vs K2 state / K2 reference, C^_EE; Gamma in c2 units; EE-own basis; Jacobi fix")
store = {}
for spec in sys.argv[1].split(','):
    J, D0 = int(spec.split(':')[0]), float(spec.split(':')[1])
    db = ED.dload(J)
    W.DSTATE = False
    g2 = Q.A.grid(J, 2, ('C2', 'C3H'))
    W.DSTATE = True
    gd = Q.A.grid(J, 2, ('C2', 'C3H'))
    W.DSTATE = False
    r2, rd = g2[D0], gd[D0]
    m = r2['m']
    sel = np.array([any(k[0] == D0 and abs(k[1] - x*m) < 1e-6*x*m for k in db) for x in r2['xs']])
    errs = {}
    for (D, nu), r in db.items():
        if D == D0:
            errs[round(nu/m, 4)] = (np.max(np.abs(r['err']))/np.max(np.abs(r['E'])), r['nmask'])
    pr(f"  J={J} Delta={D0:.0f} (m={m:.1f}, J^2/m={J*J/m:.2f}): direct state points at x = {list(r2['xs'][sel])}; a_h direct/K2 " +
       ", ".join(f"{h} {rd['a'][h]/r2['a'][h]:.6f}" for h in r2['a']) + f"; RQMC max rel err / nmask: {errs}")
    out = {}
    for e in ('++', '+-'):
        Dh = Q.A.Dsh(e, J, r2['xs'])
        ra = m*m*((rd[('C2', e)]/m**2 + Dh)*rd['a'][e]/r2['a'][e] - Dh)
        out[e] = (ra[sel], r2[('C2', e)][sel])
        # RQMC error in Gamma units: m^2 |D-envelope| x relative error (upper estimate)
        pr(f"    EE {e}: x " + " ".join(f"{x:7.3f}" for x in r2['xs'][sel]))
        pr(f"       direct: " + " ".join(f"{v:7.0f}" for v in ra[sel]))
        pr(f"       K2    : " + " ".join(f"{v:7.0f}" for v in r2[('C2', e)][sel]))
        pr(f"       m^2 |D_h| x RQMC rel err: " + " ".join(f"{m*m*abs(d)*errs.get(round(x, 4), (np.nan,))[0]:7.0f}" for x, d in zip(r2['xs'][sel], Dh[sel])))
    pr(f"    EE off: direct " + " ".join(f"{v:7.3f}" for v in rd[('C2', 'off')][sel]) + " | K2 " + " ".join(f"{v:7.3f}" for v in r2[('C2', 'off')][sel]))
    store[(J, D0)] = dict(xs=r2['xs'][sel], out=out, k2row=r2, drow=rd)
pickle.dump(store, open(os.path.join(C, 'ee1000_directcmp.pkl'), 'wb'))
open(os.path.join(C, 'ee1000_directcmp.txt'), 'w').write("\n".join(lines) + "\n")
