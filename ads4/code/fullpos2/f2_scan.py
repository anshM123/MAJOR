"""[DV] pointwise Phi2 - Phi2_hyb versus the size needed: 2 int |psi2(p)| |dPhi(Lam p)| dp  compared with N_s (T0 hyb_fast)."""
import sys, os, math, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
import numpy as np, mpmath as mp
import f2_num as F
import ls_core as L
Dg = 1000; lam = (Dg - 2)*(Dg - 1); Lam = math.sqrt(lam)
Fw = L.NFw(Dg)
def run(J, bhat, b, npts=400):
    h = '++' if b == 0 else '+-'
    m = (2*J + 1)*Lam/bhat; T = m*m + J*(J + 1)
    mfl = m/Lam
    hy = L.hyb_fast(Dg, J, mfl, h)
    ps = np.linspace(0.0005, 0.9995, npts); wq = np.full(npts, 1.0/npts)
    psi = Fw.psi2(ps)
    D = []; P2 = []; PH = []
    t0 = time.time()
    for p in ps:
        nu = Lam*p
        mp.mp.dps = int(0.75*nu + 60 + 0.1*J)
        st = F.state(J, b, T)
        a = F.Phi2(st, nu); c = F.Phi2_hyb(st, nu)
        P2.append(float(a)); PH.append(float(c)); D.append(float(a - c))
    P2, PH, D = map(np.array, (P2, PH, D))
    I_ang = 2*(psi*P2*wq).sum(); I_hyb = 2*(psi*PH*wq).sum(); Iabs = 2*(abs(psi*D)*wq).sum(); Iabs0 = 2*(abs(psi*PH)*wq).sum()
    print(f'J={J} b={b} bhat={bhat} Delta-J={1.5+math.sqrt(T+.25)-J:.0f} r={J*(J+1)/T:.4f}: N_s={hy["Ns"]:.4e} B2+cv2={hy["B2"]+hy["cv2"]:.4e}  '
          f'2intpsi*Phi2={I_ang:.4e} 2intpsi*Phihyb={I_hyb:.4e}  2int|psi dPhi|/N_s={Iabs/hy["Ns"]:.3e}  2int|psi Phihyb|/N_s={Iabs0/hy["Ns"]:.3e} '
          f'max|dPhi|={abs(D).max():.3e} at p={ps[abs(D).argmax()]:.3f}  [{time.time()-t0:.0f}s]', flush=True)
    return ps, P2, PH
if __name__ == '__main__':
    for s in sys.argv[1:]:
        J, bh, b = s.split(','); run(int(J), float(bh), int(b))
