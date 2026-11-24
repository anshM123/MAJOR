"""
gn_refdiag (GAP_NEARGAP.md §2.4): where does the K2 - K1 step of an odd-J row come from -- the state's Laplace series or the heavy
reference's (Delta_ref = 2e4)?  Per node of a K2 row (C3H and EE records, entry [-1,-1] as p4_odd_analysis.grid uses), in c2 units
(m^2 x relative): the K = 2 term of the ratio t2 = a2 + a1 e2 + a0 e4 split into the state part a2/a0 and the reference parts
(a1/a0) e2 + e4 (e2, e4 from the reference cumulants, p4_odd_analysis.pred).  Read-only on every cache.
Usage: python3 gn_refdiag.py J D[,D...]
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gn_guard  # noqa: F401
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF)
import numpy as np
os.chdir(SF)
import p4_odd_analysis as P4
J = int(sys.argv[1]); Ds = [float(v) for v in sys.argv[2].split(',')]
for kind in ('EE', 'C3H'):
    R = {}
    for d, p in ((os.path.join(SF, 'cache'), 'ee1000_p3_%s_J%d.pkl'), (os.path.join(HERE, 'cache', 'rows'), 'w1_p3_%s_J%d_gn.pkl')):
        for fn in glob.glob(os.path.join(d, p % (kind, J))):
            R.update({(k[0], k[1]): v for k, v in pickle.load(open(fn, 'rb')).items()})
    for D in Ds:
        m = P4.mphys(D, J)
        print(f"# {kind} J={J} row {D:g} (m={m:.1f}): x | state a2/a0 | ref (a1/a0)e2 + e4 | total t2/a0 | state a1/a0 | ref e2   (c2 = m^2 x)")
        for (DD, nu) in sorted(k for k in R if k[0] == D):
            if (20000.0, nu) not in R:
                continue
            st, rf = R[(D, nu)], R[(20000.0, nu)]
            a = [st[o][-1, -1] for o in ('M0', 'M2', 'M4')]; r = [rf[o][-1, -1] for o in ('M0', 'M2', 'M4')]
            k2, k4, k6 = P4.cum(r); e2, e4 = -k2, -k4 + k2*k2/2
            f = m*m
            sa2 = f*abs(a[2]/a[0]); rp = f*abs(a[1]/a[0]*e2 + e4); tt = f*abs(a[2]/a[0] + a[1]/a[0]*e2 + e4)
            x = nu/m
            if abs(x*40 - round(x*40)) < 1e-6 and round(x*40) % 4 == 0:
                print(f"  {x:.3f} | {sa2:9.0f} | {rp:9.0f} | {tt:9.0f} | {f*abs(a[1]/a[0]):9.0f} | {f*abs(e2):9.0f}")
