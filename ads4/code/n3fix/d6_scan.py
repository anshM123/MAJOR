"""n3fix/d6_scan.py: point-wise n3 (grav_on_paths) / n1 (D.evaluate, original code) ratio of the C3H (and EE) integrand at real inner
points, for the original n3 and the patched one (patch_n3c3h parts a, b, a+b), state and Delta_ref = 2e4, over x.
Also: n1 with part (a) vs n1 original (must be identical: r = +1 on n1's convention).
Usage: python3 d6_scan.py J Delta x1,x2,...  [npts]"""
import nf_guard  # noqa: F401
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc')); sys.path.insert(0, SF)
os.environ['HC_CACHE'] = os.path.join(HERE, 'cache'); os.environ.setdefault('L1_DEG', '6')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import hc_points as HC, z2_points as Z2
J = int(sys.argv[1]); Ds = float(sys.argv[2]); xs = [float(v) for v in sys.argv[3].split(',')]
npts = int(sys.argv[4]) if len(sys.argv) > 4 else 8
odd = J % 2 == 1
HC.setup(odd=odd)
import patch_n3c3h as PX
PX.enable_or(False); PX.enable_branch(False)
import n1_pw as PW, n1_direct as D, n3_flow as F, n3_hybrid as Hy
rng = np.random.default_rng(1)
Z = 0.8*rng.normal(size=(npts, 7)); Z = Z[np.linalg.norm(Z, axis=1) < 3.4]
for x in xs:
    nu = x*HC.mphys(Ds, J)
    for Dl in (Ds, 20000.0):
        c, wc, L, _ = F.saddle(nu, Dl, J); Z2.set_qsc(Dl, J, c, wc, L)
        evS = F.Evaluator(nu, Dl, J, c, wc, L, 'S', nuf=nu); F.set_branch(evS); Hy.attach_nuflow(evS)
        ct = Hy.contour(evS, Z, 3.5, keep_path=True)
        line = f"x={x:.4f} Delta={Dl:g} (Im X20(wc) sign via r test): "
        for kind in ('EE', 'C3H'):
            def n1(on):
                PX.enable_or(on)
                ev = F.Evaluator(nu, Dl, J, c, wc, L, kind, nuf=nu)
                if odd:
                    import p8_odd_pw as P8M
                    g = lambda zl: P8M.P8.build(zl, nu, Dl, J, kind, 1.0)[:3]
                else:
                    g = lambda zl: PW.P9.build_spin(zl, nu, Dl, J, ev.mons, ev.sol, ev.solt, kind)
                return D.evaluate(g, c, wc, L, Z, nu, Dl, J)[0].reshape(-1, Z.shape[0])
            ref = n1(False)
            res = {}
            for (a, b) in ((False, False), (True, False), (False, True), (True, True)):
                PX.enable_or(a); PX.enable_branch(b)
                ev = F.Evaluator(nu, Dl, J, c, wc, L, kind, nuf=nu); F.set_branch(ev)
                Fg = Hy.grav_on_paths(ev, ct, Z.shape[0])[0].reshape(-1, Z.shape[0])
                res[(a, b)] = np.median(np.abs(Fg/ref - 1))
            dn1 = np.abs(n1(True)/ref - 1).max()
            PX.enable_or(False); PX.enable_branch(False)
            line += f" {kind}: n3/n1-1 orig {res[(False, False)]:.1e} a {res[(True, False)]:.1e} b {res[(False, True)]:.1e} a+b {res[(True, True)]:.1e} | n1(a)/n1-1 {dn1:.1e};"
        print(line, flush=True)
