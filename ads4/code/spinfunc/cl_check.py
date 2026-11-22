"""
cl_check (RESULT_D_CLEANUP.md): per-node validation of this agent's direct small-nu points (cl_oopts_* 3x3 at even J, cl_pts_* 1x1 at
odd J): (a) the scalar check against the exact Mellin object (z3_nodes.scalar_check, unchanged: d_D - d_ref must vanish; Delta c2 with
jackknife error), the protocol's contour test; (b) the max relative RQMC error of the graviton elements (state, reference).
Exact-scalar records go to cache/cl_exact.pkl (initialised from hc_exact.pkl by cl_chain.py).
Usage: python3 cl_check.py J:D[:oo] [...]  -> cache/cl_check.txt (appended)
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cl_guard     # noqa: E402,F401
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import z3_nodes as Z3
C = os.path.join(HERE, 'cache')
DREF = 20000.0
Z3.EXPK = os.path.join(C, 'cl_exact.pkl')
lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


for spec in sys.argv[1:]:
    p = spec.split(':'); J, D = int(p[0]), float(p[1]); oo = len(p) > 2
    f = os.path.join(C, ('cl_oopts_EE_J%d_D%g.pkl' if oo else 'cl_pts_EE_J%d_D%g.pkl') % (J, D))
    P = pickle.load(open(f, 'rb'))
    m = Z3.mphys(D, J)
    pr(f"## {os.path.basename(f)}  (m = {m:.3f})")
    for nu in sorted(set(k[1] for k in P if k[0] == D and (DREF, k[1]) in P and k[1] < 60)):
        dd = {DD: P[(DD, nu)] for DD in (D, DREF)}
        o = Z3.scalar_check(J, nu, dd, D)
        rel = [float(np.max(np.abs(dd[DD]['err']))/np.max(np.abs(dd[DD]['E']))) for DD in (D, DREF)]
        pr(f"  nu={nu:7.3f} x={nu/m:.4f}: scalar d-d_ref Re {o['d'].real:+.4f}+-{o['err'].real:.4f} Im {o['d'].imag:+.4f}+-{o['err'].imag:.4f};"
           f" Dc2 {o['dc2']:+7.0f}+-{o['ec2']:5.0f} ({abs(o['dc2'])/max(o['ec2'], 1e-9):.1f} sigma);  graviton max rel err state {rel[0]:.3f} ref {rel[1]:.3f}")
with open(os.path.join(C, 'cl_check.txt'), 'a') as fh:
    fh.write("\n".join(lines) + "\n")
