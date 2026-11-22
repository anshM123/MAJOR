"""
cl_j48 (RESULT_D_CLEANUP.md, task 3): consistency check of the fixed-engine K2 row (1048, 48) of FOLD_WINDOW_LARGEJ §4 against DIRECT
n1 nodes (cl_points.py, CL_FWPATCH=1: the fw structure basis; state + reference 2e4 on common random numbers) at K2-grid x.
Gamma_{k,h} (c2 units, l9/q3/w2 EE conventions, fw_gamma's generalised helicity basis from the K2 order-0 records) is built twice by
fw_gamma.row (unchanged): with the direct nodes replacing the K2 nodes at the same nu, and with FW_NODIRECT=1 (pure K2).  The curvature
image is subtracted in both (sm_curvfix.curv_image), as in the chain.  Jackknife over the 8 common-seed replicas of the direct nodes.
Usage: python3 cl_j48.py  -> cache/cl_j48.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cl_guard     # noqa: E402,F401
import numpy as np
import fw_gamma as FG
import sm_curvfix as CF
C = os.path.join(HERE, 'cache')
J, D = 48, 1048.0
SKIP = [float(v) for v in os.environ.get('CL48_SKIP', '0.475').split(',') if v]   # failed reference (E00/E01 = 82 vs ~1.4)


def fwpts(kind, J_, D_=None):
    f = os.path.join(C, 'cl_pts_%s_J%d_D%g.pkl' % (kind, J_, D))
    P = pickle.load(open(f, 'rb')) if (J_ == J and os.path.exists(f)) else {}
    m = FG.Z3.mphys(D, J)
    return {k: v for k, v in P.items() if not any(abs(k[1]/m - s) < 1e-6 for s in SKIP)}


FG.fwpts = fwpts
_gb = FG.gbasis
FG.gbasis = lambda J_, D_, R: _gb(J_, D_, FG.Q._strip(dict(FG.fwk2('EE', J_)), J_)) if J_ == J else _gb(J_, D_, R)   # pure-K2 basis
out = []


def pr(s=""):
    print(s, flush=True); out.append(s)


P = fwpts('EE', J)
nus = sorted(k[1] for k in P if k[0] == D and (FG.DREF, k[1]) in P)
K2NU = sorted(set(k[1] for k in FG.fwk2('EE', J) if k[0] == D))
for nu in nus:
    assert any(abs(nu - v) < 1e-7*nu for v in K2NU), ('direct node not on the K2 grid', nu)  # basis: pure K2 order-0 records (as the chain)


def gam(jk=None, nodirect=False):
    os.environ['FW_NODIRECT'] = '1' if nodirect else '0'
    FG.CUR.clear(); FG.CUR[J] = dict(jk=jk, D=D)
    FG.Z3._SH.clear()
    m = FG.Z3.mphys(D, J)
    g = FG.Q.A.grid(J, 2, ('C2', 'C3H'), XS=[nu/m for nu in K2NU], Ds=[D])     # the full K2 x grid in both runs
    FG.CUR.clear()
    r = g[D]
    res = {}
    for kd in ('C2', 'C3H'):
        for h in ('++', '+-'):
            v = np.asarray(r[(kd, h)], float) - CF.curv_image(kd, h, J, r['m'], r['xs'], FG.Q.A.env, FG.Q.A.SIG)
            res[(kd, h)] = v
    return r['xs'], res


xs, Gd = gam()
_, Gk = gam(nodirect=True)
jk = [gam(jk=i)[1] for i in range(8)]
pr(f"# cl_j48: (1048, 48) Gamma (c2, curvature image subtracted) at the direct nodes (row = full K2 grid, a_h refitted): direct (jackknife) | K2 (fixed-engine fw_k2 row) | difference")
for i, x in enumerate(xs):
    if not any(abs(x*FG.Z3.mphys(D, J) - nu) < 1e-6*nu for nu in nus):
        continue
    s = f"  x={x:.4f} nu={x*FG.Z3.mphys(D, J):7.2f}:"
    for k in Gd:
        V = np.array([g[k][i] for g in jk]); e = np.sqrt(7/8*((V - V.mean())**2).sum())
        s += f"  {k[0]} {k[1]}: {Gd[k][i]:+7.0f}+-{e:4.0f} | {Gk[k][i]:+7.0f} | {Gd[k][i] - Gk[k][i]:+6.0f}"
    pr(s)
pr(f"# x-nodes of the direct points: {[round(nu/FG.Z3.mphys(D, J), 4) for nu in nus]}")
open(os.path.join(C, 'cl_j48.txt'), 'w').write("\n".join(out) + "\n")
