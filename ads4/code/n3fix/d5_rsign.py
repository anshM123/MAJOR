"""n3fix/d5_rsign.py: the relative sign r = n2 X40 / (n4 X20) of the two terms of o_H (l9_points.odd_structure: n2 = (X24 X20 X40^3)^-1/2,
n4 = (X24 X40 X20^3)^-1/2; on the physical sheet both equal (X24 X20 X40)^-1/2 / X40, / X20, i.e. r = +1).  Evaluated (a) with principal
logs at the crystal (TPS / Laplace convention), (b) principal at the full saddle wc (n1's convention), (c) continued from W_REF next to the
crystal to wc (n3's set_branch).  python3 d5_rsign.py J Delta x"""
import nf_guard  # noqa: F401
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc')); sys.path.insert(0, SF)
os.environ['HC_CACHE'] = os.path.join(HERE, 'cache'); os.environ.setdefault('L1_DEG', '6')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import hc_points as HC, z2_points as Z2
J = int(sys.argv[1]); x = float(sys.argv[3])
HC.setup(odd=False)
import n1_pw as PW, n1_direct as D, n3_flow as F
for Dl in (float(sys.argv[2]), 20000.0):
    nu = x*HC.mphys(float(sys.argv[2]), J)
    c, wc, L, _ = F.saddle(nu, Dl, J); Z2.set_qsc(Dl, J, c, wc, L)
    ev = F.Evaluator(nu, Dl, J, c, wc, L, 'C3H', nuf=nu); F.set_branch(ev)
    def XX(w):
        PW.set_batch(w.shape[0]); zl, lj = D.zmap([w[:, i].copy() for i in range(7)], c, pw=True)
        P = PW.A6.geom(zl, 1.0)
        return PW.S.Xfun(P['x2'], P['x0']), PW.S.Xfun(P['x4'], P['x0']), PW.S.Xfun(P['x2'], P['x4'])
    pts = {'crystal': np.zeros(7, complex), 'wc': wc}
    for name, w in pts.items():
        X20, X40, X24 = XX(w[None, :])
        n2 = (X24*X20*X40**3)**-0.5; n4 = (X24*X40*X20**3)**-0.5          # principal at that point
        print(f"Delta={Dl:g} x={x}: {name}: X20 {X20[0]:.3e} X40 {X40[0]:.3e} X24 {X24[0]:.3e}; principal r = {(n2*X40/(n4*X20))[0]:+.3f}")
    # n3 convention: saddle refs (logs 31, 32 = the two powr of odd_structure, by d4_logid)
    sr = ev.saddle_refs['grav']
    X20, X40, X24 = XX(wc[None, :])
    r3 = np.exp(-0.5*sr[31] + 0.5*sr[32])*X40/X20
    print(f"   n3 set_branch r at wc = {complex(np.ravel(r3)[0]):+.3f};  log31 check {np.abs(np.exp(sr[31]) - X24*X20*X40**3).max()/np.abs(X24*X20*X40**3).max():.1e}")
