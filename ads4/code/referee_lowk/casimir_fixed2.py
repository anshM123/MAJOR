"""
The note's Casimir check (lowk/casimir_check.py) computes 1 - z by subtraction; near rho -> 1 this loses eps/(1-z) ~ 1e-16/1e-11
relative accuracy, i.e. the quoted deep-corner residual (~5e-5 at 1-rho_b = 3e-6, identical for A and B) is the check's own roundoff.
Here: same residual with 1 - z = ((1-rho)/(1+rho))^2 computed without cancellation, for settings A, B and a window/ansatz variant.
Usage: python3 casimir_fixed.py
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'lowk')); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np
from lcblock import LCBlock
from c1_block_checks import corner_points
from casimir_check import casimir_residual
from corner import PAR_A, PAR_B

def casimir_residual_fixed(E, rs, rb, Delta, J, a, b, d=3):
    F, Fs, Fb, Fss, Fbb = E[0, 0], E[1, 0], E[0, 1], E[2, 0], E[0, 2]
    def one(rho, F1, F2):
        w = (1 + rho)/(1 - rho); thw = 2*rho/(1 - rho)**2
        z = 4*rho/(1 + rho)**2; omz = ((1 - rho)/(1 + rho))**2
        zd = w*F1; zd2 = w*w*F2 + w*thw*F1
        Dz = omz*(zd2 - zd) - (1 + a + b)*z*zd - a*b*z*F
        sc = np.maximum.reduce([np.abs(omz*zd2), np.abs(omz*zd), np.abs((1 + a + b)*z*zd)])
        return Dz, z, omz, zd, sc
    Ds, zs, omzs, zds, scs = one(rs, Fs, Fss)
    Db, zb, omzb, zdb, scb = one(rb, Fb, Fbb)
    m1 = (d - 2)/(zs - zb)*zb*omzs*zds; m2 = (d - 2)/(zs - zb)*zs*omzb*zdb
    C = Delta*(Delta - 3) + J*(J + 1)
    scale = np.maximum.reduce([scs, scb, np.abs(m1), np.abs(m2), np.abs(C/2*F)])
    return np.abs(Ds + Db + m1 - m2 - C/2*F)/scale

for (Delta, J, c1, c2) in ((44.56516128923552, 8, -0.375, 0.625),):
    pts = corner_points(Delta, J, (1.0, 3.0, 20.0), (1e-4, 1e-3, 1e-2))
    r = np.array([p[2] for p in pts]); chi = np.array([p[3] for p in pts])
    xs, xb = r*np.exp(-chi), r*np.exp(chi)
    print(f"# (Delta,J)=({Delta},{J}) (c1,c2)=({c1},{c2})")
    big = LCBlock(Delta, J, c1, c2, **PAR_A)
    res = {}
    for name, par in (("A", PAR_A), ("B", PAR_B), ("M12", dict(PAR_A, M=12, Malt=5)), ("Q3000", dict(PAR_A, Q=3000)), ("win_lo", dict(PAR_A, qlo=1000))):
        lb = big if name == "A" else LCBlock(Delta, J, c1, c2, coeffs=big.C, **par)
        E = lb.E(r, chi)
        res[name] = (casimir_residual(E, xs, xb, Delta, J, c1, c2), casimir_residual_fixed(E, xs, xb, Delta, J, c1, c2), lb.series(r, chi))
    for p, (t, u, _, _) in enumerate(pts):
        dd = {n: max(abs(res[n][2][p, k, l]/res['A'][2][p, k, l] - 1) for k in range(5) for l in range(5 - k)) for n in res if n != 'A'}
        print(f"  t={t:4.1f} u={u:.0e} 1-rho_b={1-xb[p]:.1e}: Casimir (note's formula) A {res['A'][0][p]:.1e} | fixed 1-z: " +
              " ".join(f"{n} {res[n][1][p]:.1e}" for n in res) + " | max|S_X/S_A-1|: " + " ".join(f"{n} {v:.1e}" for n, v in dd.items()), flush=True)
