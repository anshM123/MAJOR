"""n3fix/d1_pointwise.py: point-wise comparison of the graviton integrand (kinds EE, C3H) on the n3 contour (grav_on_paths: logs replayed
along the stored ray/flow paths) and through n1's evaluator (D.evaluate: two-level ray continuation from the saddle), at REAL zeta inside
the inner disk, where the two contours coincide (xi = zeta), with the same saddle (c, wc, L).
Usage: python3 d1_pointwise.py J Delta x [npts] [scale]"""
import nf_guard  # noqa: F401
import os, sys, time
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
npts = int(sys.argv[4]) if len(sys.argv) > 4 else 24
sc = float(sys.argv[5]) if len(sys.argv) > 5 else 0.8
odd = J % 2 == 1
HC.setup(odd=odd)
import n1_pw as PW, n1_direct as D, n3_flow as F, n3_hybrid as Hy, l5_spin as S
nu = x*HC.mphys(Dl, J)
c, wc, L, _ = F.saddle(nu, Dl, J)
Z2.set_qsc(Dl, J, c, wc, L)
evS = F.Evaluator(nu, Dl, J, c, wc, L, 'S', nuf=nu); F.set_branch(evS); Hy.attach_nuflow(evS)
rng = np.random.default_rng(1)
Z = sc*rng.normal(size=(npts, 7))
Z = Z[np.linalg.norm(Z, axis=1) < 3.4]
ct = Hy.contour(evS, Z, 3.5, keep_path=True)
print(f"J={J} Delta={Dl} x={x} nu={nu:.4f}; {Z.shape[0]} inner points; xi==zeta: {np.abs(ct['xi'] - Z).max():.1e}", flush=True)
for kind in ('EE', 'C3H'):
    if odd:
        ev = F.Evaluator(nu, Dl, J, c, wc, L, kind, nuf=nu)
    else:
        ev = F.Evaluator(nu, Dl, J, c, wc, L, kind, nuf=nu)
    F.set_branch(ev)
    Fg, Pgc, pref0, okg = Hy.grav_on_paths(ev, ct, Z.shape[0])
    if odd:
        import p8_odd_pw as P8M
        gfun = lambda zl: P8M.P8.build(zl, nu, Dl, J, kind, 1.0)[:3]
    else:
        gfun = lambda zl: PW.P9.build_spin(zl, nu, Dl, J, ev.mons, ev.sol, ev.solt, kind)
    val, ok, _ = D.evaluate(gfun, c, wc, L, Z, nu, Dl, J)
    a = Fg.reshape(-1, Z.shape[0]); b = val.reshape(-1, Z.shape[0])
    r = a/b
    print(f"{kind}: pref0 {np.round(pref0.ravel()[:4], 6)}; n3/n1 ratio per point (element 0): ", flush=True)
    for i in range(min(Z.shape[0], 12)):
        print(f"   |z|={np.linalg.norm(Z[i]):.2f}: n3 {a[0, i]:+.3e}  n1 {b[0, i]:+.3e}  ratio {r[0, i]:+.4f}  max|ratio-r0| over elems "
              f"{np.abs(r[:, i] - r[0, i]).max():.1e}", flush=True)
    print(f"  {kind}: median |ratio - 1| {np.median(np.abs(r - 1)):.2e}; max {np.abs(r - 1).max():.2e}; ok n3 {okg.sum()} n1 {ok.sum()}", flush=True)
