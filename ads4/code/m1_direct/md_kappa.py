"""md_kappa.py -- J = 0: fit the direct D' image to  c [z2/3 + kappa 2 y2/T]  (kappa = 1 is the ladder; kappa = 0 leading Regge)
and to c [z2/3 + 2y2/T + mu w(nu)/T^2] style deformations; reference state Delta_ref with the ladder P (grid-matched).
Usage: python3 md_kappa.py Dref numax"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import md_compare_J0 as MC


def main(Dref, numax):
    import glob
    Ds = sorted(set(float(os.path.basename(f).split('_')[1][1:]) for f in glob.glob(os.path.join(HERE, 'out', 'dq_D*_E2_*.npz'))))
    RR = {g: (n, I) for g, n, I in MC.best(Dref)}
    for D in Ds:
        if D == Dref:
            continue
        for g, nu, I in MC.best(D):
            if g not in RR:
                continue
            nuR, IR = RR[g]
            m = np.isin(np.round(nu, 6), np.round(nuR, 6)) & (nu <= numax)
            nu1 = nu[m]; I1 = I[m]; I2 = np.array([IR[np.argmin(np.abs(nuR - v))] for v in nu1])
            s = nu1*nu1 + 0.25; T = (D - 1)*(D - 2)
            P2 = MC.P(nu1, Dref)
            a = I1*P2
            B = np.stack([I2*MC.z2(s)/3, I2*2*MC.y2(s)/T], axis=1)
            x, *_ = np.linalg.lstsq(B, a, rcond=None)
            res = a - B @ x
            env = np.sqrt(np.mean(np.abs(a)**2))
            print(f'  Delta={D:g} grid {g}: kappa = {(x[1]/x[0]).real:+.5f}{(x[1]/x[0]).imag:+.1e}i  max|res|/rms = {np.max(np.abs(res))/env:.1e}'
                  f'  (fixed kappa=1: {np.max(np.abs(MC.fit(nu1, I1, I2, MC.P(nu1, D), P2, np.ones_like(nu1))[1])):.1e})')


if __name__ == '__main__':
    main(float(sys.argv[1]), float(sys.argv[2]))
