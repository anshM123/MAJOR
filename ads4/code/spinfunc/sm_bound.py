"""
sm_bound (SMOOTH_BOUND.md §1): weights of the possible non-oscillating (smooth) beyond-Regge structures in the psi-integrated
Omega^EE action of a near-gap state, and the resulting tolerances, on the hybrid states of ee1000_common.hyb (k = 2 part, as
SUBLEADING_REGGE's tolerance table; p-rule 512 nodes).  Per state (gap mass, Delta = Delta_g + J) and helicity:
  W_all     = sum p2 env2 / m^2 / N_s                (uniform smooth Gamma = 1 c2)  -> tol_all = 0.9 lam/W_all  (reproduces SR §0.1)
  W_F(nuF)  = same restricted to nu < nuF            (fold region, nu* ~ 11.3; nuF = 20, 30)
  W_F^an    = X0 nuF^2 /(Lam^2 m^2 N_s) (env(0) = 2, psi2 ~ X0 p, X0 = 188.9: analytic small-p form)
  W_E       = sum p2 env2 f_E / m^2 / N_s,  f_E = ((1 - x_max^2)/(1 - x^2))^3  (edge-shaped smooth term, = 1 at x_max)
  E_osc     = endpoint estimate of an OSCILLATING edge term A (1-x_max^2)^3/(1-x^2)^3 cos(phi): 5! |C_psi| env(x_max) A /(m^2 phi'^6 N_s),
              phi' = bhat/sqrt(1 - x_max^2) per unit p, C_psi = lim psi2/(1-p)^5 = X(1)
Also lam0 = base/N_s.  Usage: python3 sm_bound.py Dg J1,J2,..  [npn]  -> stdout (cache/sm_bound_<Dg>.txt via tee)
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import ee1000_common as E

if __name__ == '__main__':
    Dg = int(sys.argv[1]); Js = [int(v) for v in sys.argv[2].split(',')]
    npn = int(sys.argv[3]) if len(sys.argv) > 3 else 512
    F = E.Q.NearGap(Dg).F
    X0 = float(np.polynomial.polynomial.polyval(0.0, F.X)); X1 = float(np.polynomial.polynomial.polyval(1.0, F.X))
    Lam = E.lam_of(Dg)
    print(f"# sm_bound Dg={Dg}: psi2 = (1-p)^5 p X(p), X(0) = {X0:.2f}, X(1) = {X1:.2f}; Lam = {Lam:.1f}")
    for J in Js:
        mfl = E.mgap(Dg, J)
        out = []
        for h in ('++', '+-'):
            st = E.hyb(Dg, J, mfl, h, npn)
            x, m, Ns = st['x'], st['m'], st['Ns']; nu = x*m
            env = E.Q.A.env(h, x, 2)
            w = st['p2']*env/m**2/Ns
            lam0 = st['base']/Ns
            xm = 1/mfl
            fE = ((1 - xm**2)/(1 - x**2))**3
            WF = {nf: w[nu < nf].sum() for nf in (20, 30, 50)}
            WFan = {nf: X0*nf**2/(Lam**2*m**2*Ns) for nf in (20, 30)}
            WE = (w*fE).sum()
            bh = (2*J + 1)/mfl
            php = bh/np.sqrt(1 - xm**2)
            Eosc = 120*abs(X1)*float(E.Q.A.env(h, np.array([xm]), 2)[0])/(m**2*php**6*Ns)
            out.append((h, lam0, w.sum(), WF, WFan, WE, Eosc, Ns))
        bh = (2*J + 1)/mfl
        print(f"J={J:3d} m/L={mfl:.4f} m={mfl*Lam:.1f} bhat={bh:6.1f} beta={2*np.arcsinh((J+.5)/(mfl*Lam)):.4f} x_max={1/mfl:.3f}")
        for h, lam0, Wa, WF, WFan, WE, Eosc, Ns in out:
            print(f"   {h}: lam0 {lam0:+.3f} N_s {Ns:.2e} | W_all {Wa:.2e} tol_all {0.9*lam0/Wa:8.2f} c2 | fold W_F(20) {WF[20]:.2e} "
                  f"(an {WFan[20]:.2e}) tol {0.9*lam0/WF[20]:9.0f}; W_F(30) {WF[30]:.2e} (an {WFan[30]:.2e}) tol {0.9*lam0/WF[30]:8.0f}; "
                  f"W_F(50) {WF[50]:.2e} tol {0.9*lam0/WF[50]:7.0f} | edge W_E {WE:.2e} tol A_E {0.9*lam0/WE:9.1f} c2; osc-edge per unit A {Eosc:.1e}")
