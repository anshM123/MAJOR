"""
c1: calibration of the position-space (real-section) moment quadrature of e2e/moments.py at the near-gap states, on the SCALAR skeleton,
where the exact answer is known (exactskel/skelexact.py, Mellin/Barnes, [DV] to 1e-10).
   R_pos(k,nu) = Pihat_{k,nu}[skeleton block] / leading prediction (e2e/t3_scalar_mixed.prediction)   (same grid, N, mask as gwaves)
   R_ex(k,nu)  = m^{2k} C_{k,nu}(exact) / (2 P_{-1/2+i nu}(cosh beta))
Both are 'exact / leading Regge formula'; R_pos - R_ex measures quadrature + cancellation error of the position-space route vs nu.
Usage: python3 c1_calib.py Dg mfl J Df   -> cache/c1_<Dg>_<mfl>_<J>_Df<Df>.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ('../e2e', '../exactskel', '../omegau'):
    sys.path.insert(0, os.path.join(HERE, p))
import numpy as np, mpmath as mp
from moments import ScalarMoments
from t3_scalar_mixed import prediction
from skelexact import ExactState
from gwaves import state
import remk
mp.mp.dps = 25
Dg, mfl, J, Df = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
Lam, m, Dv = state(Dg, mfl, J)
from common import gl01
pn, pw = gl01(64)
NUS = [0.0, 1.0, 3.0] + [float(Lam*p) for p in pn if Lam*p < 25]
al, alb = Dv - J - 2, Dv + J - 1
eta = (al*al + alb*alb)/(2*al*alb)
es = ExactState(Dv, J, Df); ex = es.values(np.array(NUS))
Pd = remk.Pders_vec(np.array(NUS), eta, nmax=0)[0]
mom = ScalarMoments(Dv, J, Df, Df, 8000)
out = [f"# c1 calibration: Dg={Dg} m/Lam={mfl} J={J} Delta={Dv:.6f} m={m:.4f} beta={np.arccosh(eta):.5f} Df={Df}",
       "#  k    nu      x=nu/m      R_ex            R_pos           R_pos-R_ex"]
for k in (2, 3, 4, 5, 6, 7):
    for i, nu in enumerate(NUS):
        sg, lg = mom.moment(k, nu)
        Rp = float(sg*mp.exp(mp.mpf(lg))/prediction(Dv, J, Df, Df, k, nu).real)
        Re = float(ex[k][i]/(2*Pd[i]))
        out.append(f"  {k}  {nu:8.4f}  {nu/m:8.5f}  {Re:+.10f}  {Rp:+.10f}  {Rp-Re:+.3e}")
txt = "\n".join(out); print(txt)
open(os.path.join(HERE, 'cache', 'c1_%d_%.3f_%d_Df%g.txt' % (Dg, mfl, J, Df)), 'w').write(txt + "\n")
