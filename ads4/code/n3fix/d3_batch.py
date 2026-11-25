"""n3fix/d3_batch.py: is the point-wise C3H prefactor really point-wise?  build_spin on the batch [saddle, p1..pn] vs [saddle, p_i] alone
(plain mode: principal logs relative to index 0), and the intermediate odd structure / q.  Usage: python3 d3_batch.py J Delta x"""
import nf_guard  # noqa: F401
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc')); sys.path.insert(0, SF)
os.environ['HC_CACHE'] = os.path.join(HERE, 'cache')
os.environ.setdefault('L1_DEG', '6')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import hc_points as HC
import z2_points as Z2

J = int(sys.argv[1]); Dl = float(sys.argv[2]); x = float(sys.argv[3])
HC.setup(odd=False)
import n1_pw as PW, n1_direct as D, n3_flow as F, l5_spin as S5
nu = x*HC.mphys(Dl, J)
c, wc, L, _ = F.saddle(nu, Dl, J)
Z2.set_qsc(Dl, J, c, wc, L)
mons, sol, _ = S5.conserved(Dl, J); _, solt, _ = S5.conserved(3.0 - Dl, J)
rng = np.random.default_rng(1)
Z = 0.8*rng.normal(size=(4, 7))


def pref_at(Zs, kind):
    Wp = wc[None, :] + np.vstack([np.zeros((1, 7)), Zs]) @ L.T
    PW.set_batch(Wp.shape[0]); zl, lj = D.zmap([Wp[:, i].copy() for i in range(7)], c, pw=True)
    lamF, G, pref = PW.P9.build_spin(zl, nu, Dl, J, mons, sol, solt, kind)
    return pref[..., 1:]


for kind in ('C2', 'C3H'):
    full = pref_at(Z, kind)
    for i in range(Z.shape[0]):
        one = pref_at(Z[i:i + 1], kind)
        print(f"{kind} point {i}: batch-of-{Z.shape[0]} vs alone: max rel diff {np.abs(full[..., i] - one[..., 0]).max()/np.abs(one).max():.2e}")
