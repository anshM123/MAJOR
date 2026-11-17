"""c2: residual flip interference of the superconvergent even moments C^SC_4, C^SC_6 ((1,3) j'=0 celestial, (2,4) S_ff flip-free) inside Omega^SC (M02),
relative to N_s, on the near-gap states of c1; and of the forward piece d_{nu^2} C^SC_4|_0.  Only one-pair flips at (1,3) survive (flip13 x pres24).
Bound per unit |c|^2: |rho13(nu)| |K_{4,0}(beta)| / 2  (|c_s cbar_o| <= (|c_s|^2+|c_o|^2)/2; O(1) bookkeeping factors not fixed) [M]."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'omegau')); sys.path.insert(0, HERE)
import numpy as np
from common import load, state, EPS
from rcommon import peval, flat_parts
from c1_flip_contamination import K, pn, pw
def rho_sc(nu, k):
    c = {4: (210, 216, 263), 6: (1890, 536, 2023)}[k]
    return c[0]*(2*nu - 5j)*(2*nu - 1j)/((2*nu + 3j)*(2*nu + 7j)*(16*nu**4 + c[1]*nu**2 + c[2]))
F, PH = load('M02')
print("c2: Omega^SC (M02): residual one-pair flip of C^SC_4, C^SC_6 / N_s  and forward d_{nu^2}C^SC_4|_0 flip / N_s")
for Dg in (100, 1000):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    for h in ('++', '+-'):
        for mfl in (1.2, 1.4, 1.6, 2.0):
            u = 1/mfl**2
            Js, B2, B3, FW = flat_parts(F, h, u, 12)
            for J in ((0, 2, 4, 6) if h == '++' else (4, 6, 8)):
                st = state(Dg, Lam, mfl, J)
                if st is None: continue
                b, m = st['beta'], st['m']
                Ns = abs(B2[J - Js[0]]) + abs(B3[J - Js[0]]) + abs(FW[J - Js[0]])
                fl = 0.0
                for i, p in enumerate(pn):
                    nu = Lam*p; k40 = abs(K(nu, b, 4, 0))
                    for k in (4, 6):
                        fl += pw[i]*u**(k - 2)*abs(peval(PH[k], p))*2*EPS[k]*abs(rho_sc(nu, k))*k40/2
                # forward: -a u^3 d_y [eps4 m^8 C4] = a u^3 m^2 d_{nu^2}[...]; flip part 2 rho13 F_4 X, bound with X <= 1/2
                hstep = 0.05
                g = lambda nu: (rho_sc(nu, 4)*K(nu, b, 4, 0))
                d = (g(np.sqrt(2*hstep)) - g(np.sqrt(hstep)))/hstep      # d/d(nu^2), forward difference away from 0
                fw = abs(F.a)*u**3*m*m*abs(d)
                print(f"  {Dg:5d} {h} {mfl:4.2f} J={J:2d}  smeared flip/N_s = {fl/Ns:.3g}   forward flip/N_s = {fw/Ns:.3g}", flush=True)
