"""edger/er_diag.py (EDGE_R.md §2): why the n1 re-centred contour fails at nu/Delta >~ 0.94.  Scalar crystal only (no graviton), on a
state row: saddle, Hessian spectrum, and the fraction / kind of masked proposal points (wall: Re(Phi - Phi_0) > 2; winding > WMAX).
Usage: python3 er_diag.py J Delta x1,x2,...   -> stdout"""
import er_guard  # noqa
import os, sys
sys.dont_write_bytecode = True
SF = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'spinfunc'); sys.path.insert(0, SF)
os.environ.setdefault('L1_DEG', '6')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import n1_direct as D, n1_pw as PW, l1_tps as T
import z2_points as Z2, fw_frame as FW


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


if __name__ == '__main__':
    J = int(sys.argv[1]); Dl = float(sys.argv[2]); m = mphys(Dl, J)
    FW.patch_qrec_pw(PW, SF); Z2.patch()
    rng = np.random.default_rng(1)
    for x in [float(v) for v in sys.argv[3].split(',')]:
        nu = x*m
        c, wc, L = D.full_saddle(nu, Dl, J)
        Z2.set_qsc(Dl, J, c, wc, L)
        wv = [T.var(i, wc[i]) for i in range(T.NV)]
        g, H = T.grad_hess(D.phi_w(wv, c, nu, Dl, J))
        ev = np.linalg.eigvals(-H)
        xi = rng.normal(size=(4096, 7))*1.2
        val, ok, dph = D.evaluate(None, c, wc, L, xi, nu, Dl, J)
        wind = PW.WINDP[1:]
        print(f"x={x:.4f} nu/Delta={nu/Dl:.4f}: |Hess eig| min/max {np.abs(ev).min():.3g}/{np.abs(ev).max():.3g}; "
              f"masked {np.mean(~ok):.3f} (wall {np.mean(dph >= 2):.3f}, wind {np.mean(wind >= D.WMAX):.3f}); c = {np.round(c, 4)}; "
              f"wc = {np.round(wc, 3)}", flush=True)
