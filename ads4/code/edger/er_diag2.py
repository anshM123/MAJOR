"""edger/er_diag2.py: scalar crystal on the n1 contour at a near-edge node: plain Gaussian proposal (scale s) vs the adaptive
moment-matched one (hc_points' default); relative RQMC error and masked count.  Usage: python3 er_diag2.py J Delta x mq"""
import er_guard  # noqa
import os, sys, time
sys.dont_write_bytecode = True
SF = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'spinfunc'); sys.path.insert(0, SF)
os.environ.setdefault('L1_DEG', '6')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import n1_direct as D, n1_pw as PW
import z2_points as Z2, fw_frame as FW
from er_diag import mphys

if __name__ == '__main__':
    J = int(sys.argv[1]); Dl = float(sys.argv[2]); x = float(sys.argv[3]); mq = int(sys.argv[4])
    FW.patch_qrec_pw(PW, SF); Z2.patch()
    m = mphys(Dl, J); nu = x*m
    c, wc, L = D.full_saddle(nu, Dl, J)
    Z2.set_qsc(Dl, J, c, wc, L)
    for s in (0.8, 1.0, 1.2):
        t = time.time()
        E, err, nm = D.integrate_qmc(None, c, wc, L, nu, Dl, J, mq, 8, s)
        print(f"x={x} plain s={s}: E={complex(np.ravel(E)[0]):.6e} rel err {float(np.max(np.abs(err)/np.abs(E))):.2e} nmask {nm} [{time.time()-t:.0f}s]", flush=True)
    t = time.time()
    mu, Cv = D.adapt(None, c, wc, L, nu, Dl, J, n=8192)
    E, err, nm = D.integrate_ad(None, c, wc, L, nu, Dl, J, mu, Cv, m=mq, R=8)
    print(f"x={x} adaptive: E={complex(np.ravel(E)[0]):.6e} rel err {float(np.max(np.abs(err)/np.abs(E))):.2e} nmask {nm}; |mu| {np.abs(mu).max():.2f}, "
          f"eig C {np.round(np.linalg.eigvalsh(Cv), 2)} [{time.time()-t:.0f}s]", flush=True)
