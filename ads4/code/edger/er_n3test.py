"""edger/er_n3test.py (EDGE_R.md §2): the n3 hybrid wall-free contour (Phi-thimble beyond rho1) for the scalar crystal at near-edge
nodes, with the QSC / Q_J-recurrence patches of hc_points.  Usage: python3 er_n3test.py J Delta x[,x..] mq rho1   -> stdout"""
import er_guard  # noqa
import os, sys, time
sys.dont_write_bytecode = True
SF = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'spinfunc'); sys.path.insert(0, SF)
os.environ.setdefault('L1_DEG', '6')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import n1_pw as PW, z2_points as Z2, fw_frame as FW
import n3_hybrid as NH, n3_flow as F
from er_diag import mphys

if __name__ == '__main__':
    J = int(sys.argv[1]); Dl = float(sys.argv[2]); xs = [float(v) for v in sys.argv[3].split(',')]
    mq = int(sys.argv[4]); rho1 = float(sys.argv[5])
    FW.patch_qrec_pw(PW, SF); Z2.patch()
    m = mphys(Dl, J)
    for x in xs:
        nu = x*m; t = time.time()
        c, wc, L, _ = F.saddle(nu, Dl, J)
        Z2.set_qsc(Dl, J, c, wc, L)
        for ad in (False, True):
            try:
                r = NH.integrate(nu, Dl, J, 'S', rho1=rho1, m=mq, R=8, adaptive=ad)
                print(f"x={x} nu/D={nu/Dl:.4f} n3 rho1={rho1} adaptive={ad}: Es={r['Es']:.6e} rel err {abs(r['errs']/r['Es']):.2e} "
                      f"nbad {r['diag']['nbad']}/{r['diag']['n']} [{time.time()-t:.0f}s]", flush=True)
            except Exception as e:
                print(f"x={x} adaptive={ad}: failed {e!r} [{time.time()-t:.0f}s]", flush=True)
