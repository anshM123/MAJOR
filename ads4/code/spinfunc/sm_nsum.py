"""
sm_nsum (SMOOTH_BOUND.md §2.2): n-sum convergence of the Arb Laplace/Barnes exact scalar (h13_laplace.Laplace3H) at near-gap states
of Delta_g = 300 with large J (J/tau >~ 0.13), where h13's NCONF spread is large (sm_exact: (364, 64) n-sum err ~ 40 N_s).
m^4 C_2 at nu/Lam = 0.1, 0.3, 0.6, 0.85, 0.95 for several (N0, nq, L-factor, prec, EM order).
Usage: python3 sm_nsum.py Delta J [Dg]   -> stdout (cache/sm_nsum.txt via tee)
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.join(HERE, '..', 'exactskel'); sys.path.insert(0, EX)
import numpy as np
import h13_exact as H
import h13_laplace

D, J = float(sys.argv[1]), int(sys.argv[2]); Dg = int(sys.argv[3]) if len(sys.argv) > 3 else 300
Df = 3.0; Lam = H.Lam_of(Dg)
nus = Lam*np.array([0.1, 0.3, 0.6, 0.85, 0.95])
H.mack.mack_P_fast(D, J)
print(f"# sm_nsum Delta={D} J={J} Dg={Dg}: m^4 C_2 at nu/Lam = 0.1 0.3 0.6 0.85 0.95", flush=True)
MF = float(os.environ.get('SM_MFAC', '1'))       # multiplies the w-series length M (h13 default M = 3 nu^2/(2 Delta) + 6 sqrt(nu+1) + 60)
CONFS = [(200, 960, 1.0, 400, 1), (320, 960, 0.5, 400, 1), (120, 960, 1.0, 400, 1), (80, 960, 1.0, 400, 1), (200, 960, 1.0, 800, 1),
         (200, 960, 1.0, 400, 3), (600, 960, 1.0, 800, 3), (200, 1920, 2.0, 400, 3)]
if os.environ.get('SM_CONFS'):          # "N0:nq:Lf:prec:em,..."
    CONFS = [tuple(float(v) if i == 2 else int(v) for i, v in enumerate(c.split(':'))) for c in os.environ['SM_CONFS'].split(',')]
for N0, nq, Lf, prec, em in CONFS:
    h13_laplace.EM_ORDER = em
    t0 = time.time()
    Ld = max(4.0*N0, 0.5*D**2/max(2*Df - 1.5, 0.5))*Lf
    numax = float(nus.max())
    M = int(MF*(3*numax**2/(2*D) + 6*np.sqrt(numax + 1) + 60))
    L3 = h13_laplace.Laplace3H(D, J, Df, ks=(2,), prec=prec, N0=N0, nq=nq, L=Ld, numax=numax, M=M)
    v = [L3.value_fast(2, nu) for nu in nus]
    print(f"  M={M} N0={N0:4d} nq={nq} Lf={Lf} prec={prec} EM={em}: " + ' '.join('%+.10e' % float(a.mid()) for a in v) +
          f"  [rad {max(float(a.rad())/abs(float(a.mid())) for a in v):.0e}, {time.time()-t0:.0f}s]", flush=True)
