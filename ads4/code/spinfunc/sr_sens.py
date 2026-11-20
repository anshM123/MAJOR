"""
sr_sens (SUBLEADING_REGGE.md §3): integrated sensitivity of the per-state Omega^EE action to beyond-leading-Regge corrections of a
given x-structure, on the hybrid states of ee1000_common.hyb (flat NF02 + leading-Regge curvature).  For a correction
Gamma(x) = (2J+1)^2 * B(x)   (c2 units; B a basis function) the diagonal change is dF = sum p2 env_h Gamma/m^2 (C2 part).
Basis: oscillating  c_k:  x^k D_h(x),  s_k: x^k * x D_h'(x)/(2J+1),  smooth r_k: x^k   (k = 0, 2),  D_h the LO shape.
Output: dF/N_s per unit B, per state; and the base lam0/N_s.   Usage: python3 sr_sens.py Dg J1,J2,..  [npn]
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import ee1000_common as E


def basis(h, J, x):
    D = E.Dc_fixed(h, J, x); hh = 1e-7
    Dp = (E.Dc_fixed(h, J, x + hh) - E.Dc_fixed(h, J, x - hh))/(2*hh)*x/(2*J + 1)
    return {'c0': D, 'c2': x*x*D, 's0': Dp, 's2': x*x*Dp, 'r0': np.ones_like(x), 'r2': x*x}


if __name__ == '__main__':
    Dg = int(sys.argv[1]); Js = [int(v) for v in sys.argv[2].split(',')]
    npn = int(sys.argv[3]) if len(sys.argv) > 3 else 512
    for J in Js:
        mg = E.mgap(Dg, J)
        for mfl in [mg] + [v for v in (1.1, 1.2, 1.3, 1.5, 2.0) if v > mg]:
            out = []
            for h in ('++', '+-'):
                st = E.hyb(Dg, J, mfl, h, npn)
                x, m = st['x'], st['m']
                B = basis(h, J, x)
                env = Q_env = E.Q.A.env(h, x, 2)
                sens = {k: (st['p2']*env*(2*J + 1)**2*v/m**2).sum()/st['Ns'] for k, v in B.items()}
                out.append((h, st['base']/st['Ns'], st['Ns'], sens))
            bh = 2*J + 1
            print(f"J={J:3d} m/L={mfl:.3f} bhat={bh/mfl:6.1f} " + " | ".join(
                f"{h} lam0 {l0:+.3f} Ns {Ns:.2e} " + " ".join(f"{k}:{v:+.2e}" for k, v in s.items()) for h, l0, Ns, s in out), flush=True)
