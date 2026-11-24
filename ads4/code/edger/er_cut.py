"""edger/er_cut.py (EDGE_R.md §2): the scalar crystal on the n1 contour with a radial cut |xi| <= rc (plain Gaussian proposal, scale s),
vs its Laplace series S0 + S2 + S4 at the main saddle (n1_direct.scalar convention), as a function of rc: a plateau in rc means the
hump (P_J growth at complex angle, er_diag) is separated from the main-saddle contribution by a deep valley.
Usage: python3 er_cut.py J Delta x[,x..] mq s rc1,rc2,..   -> stdout"""
import er_guard  # noqa
import os, sys, time
sys.dont_write_bytecode = True
SF = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'spinfunc'); sys.path.insert(0, SF)
os.environ.setdefault('L1_DEG', '6')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import n1_direct as D, n1_pw as PW, l1_tps as T, l1_crystal as Cr
import z2_points as Z2, fw_frame as FW
from er_diag import mphys


def laplace(nu, Dl, J, c, wc, L):
    """Laplace series of the same integral (w-contour, scalar pref form) at the full saddle: (S0, S2, S4) / value at wc."""
    wv = [T.var(i, wc[i]) for i in range(T.NV)]
    zv, lj = D.zmap(wv, c)
    lamF, G, Q = Cr.build(zv, nu, Dl, J, pref=True)
    Ph = lamF + G + lj
    # move the (J/2) log Z^2W^2 into the exponent consistently with n1: pref Q_J contains |Z|^J|W|^J -> keep as is (pref form)
    Lw = T.linsub(wc, L)
    zv2, lj2 = D.zmap(Lw, c)
    lamF2, G2, Q2 = Cr.build(zv2, nu, Dl, J, pref=True)
    S = np.array(T.expand_pref(lamF2, G2 + lj2, Q2))
    return S/Q2[0]


if __name__ == '__main__':
    J = int(sys.argv[1]); Dl = float(sys.argv[2]); xs = [float(v) for v in sys.argv[3].split(',')]
    mq = int(sys.argv[4]); s = float(sys.argv[5]); rcs = [float(v) for v in sys.argv[6].split(',')]
    FW.patch_qrec_pw(PW, SF); Z2.patch()
    m = mphys(Dl, J)
    for x in xs:
        nu = x*m
        c, wc, L = D.full_saddle(nu, Dl, J)
        Z2.set_qsc(Dl, J, c, wc, L)
        PW.set_batch(1)
        try:
            lap = laplace(nu, Dl, J, c, wc, L)
            ls = f"Laplace/pref0: S0 {lap[0]:.4f} S0+S2 {lap[0]+lap[1]:.4f} S0+S2+S4 {lap[:3].sum():.4f}"
        except Exception as e:
            ls = f"Laplace failed ({e})"
        zl, lj = D.zmap([np.array([wc[i]]) for i in range(7)], c, pw=True)
        q0 = PW.Cr.build(zl, nu, Dl, J, pref=True)[2][0]
        print(f"x={x} nu/Delta={nu/Dl:.4f}: {ls}", flush=True)
        for rc in rcs:
            t = time.time()
            E, err, nm = D.integrate_qmc(None, c, wc, L, nu, Dl, J, mq, 8, s, cut=rc if rc > 0 else None)
            E = complex(np.ravel(E)[0])/complex(q0)*(2*np.pi)**0 ; e = float(np.ravel(err)[0])/abs(complex(q0))
            print(f"   cut {rc:4.1f}: E/pref0 = {E.real:+.5f}{E.imag:+.5f}i +- {e:.1e} (nmask {nm}) [{time.time()-t:.0f}s]", flush=True)
