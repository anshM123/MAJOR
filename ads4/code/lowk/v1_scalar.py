"""
Validation on scalar partial waves (no spinning operators), <phi_1 phi_2 phi_1 phi_2>, s-channel (12) -> (34), k = 2, 3.
  exact:   (-1)^J Pihat_{k,nu}[u'^{-(D1+D2)/2} g^{a,b}]   on the (u, t) grid of ads4/code/e2e/moments.py, with
           (i) the radial series alone (as in e2e/t3_scalar_mixed.py: all grid points, N levels, no convergence mask),
           (ii) the radial series on the converged points + lightcone/tail representation (lcblock.LCBlock) on the others;
  predicted: the derived scalar heavy formula (SPINNING_E2E.md sect. 4.1; e2e/t3_scalar_mixed.prediction), nothing fitted.
Cases:  (D1, D2) = (3, 1.75): corner weight u^{2 D2 + k - 3/2} = u^{4} (k = 2): the radial series alone is converged
        (consistency check of (ii) against the published t3 numbers);
        (D1, D2) = (1.75, 0.5): corner weight u^{3/2} (k = 2), u^{5/2} (k = 3), i.e. exactly the corner weight of the spinning
        k = 2, 3 moments: the radial series alone is not converged; (ii) must give m^2 (R - 1) -> const.
Usage: python3 v1_scalar.py [quick]       Output: out_v1_scalar.txt
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np, mpmath as mp
from moments import make_grid, frame_data, series_converged
from formulas import Fnu_np
from blocks import CBlock
from lcblock import LCBlock
from t3_scalar_mixed import prediction
from cases import cases_fixed_beta, info
mp.mp.dps = 25

class ScalarLC:
    def __init__(self, Delta, J, D1, D2, N, lcpar=None, use_lc=True, itail_max=1e-13):
        self.Delta, self.J = Delta, J
        g = make_grid(Delta, J); self.g = g
        ok = g['ok']
        fd = frame_data(g['U'][ok], g['S'][ok])
        b = (D1 - D2)/2
        assert b > 0
        blk = CBlock(Delta, J, b, b, N)
        conv = series_converged(fd['rs'], fd['chi'], J, N) if use_lc else np.ones(fd['rs'].size, bool)
        Sr = np.zeros(fd['rs'].size)
        Sr[conv] = blk.series(fd['rs'][conv], fd['chi'][conv], sig=+1)[:, 0, 0]
        del blk
        self.nlc = 0; self.nleft = 0
        if use_lc and (~conv).any():
            lcpar = dict(N=10000, I=1500, Q=8000, M=10) if lcpar is None else lcpar
            lb = LCBlock(Delta, J, b, b, **lcpar)
            idx = np.where(~conv)[0]
            S = lb.series(fd['rs'][idx], fd['chi'][idx], kmax=0, lmax=0)[:, 0, 0]
            good = lb.i_tail < itail_max
            Sr[idx[good]] = S[good]; conv[idx[good]] = True
            self.nlc = int(good.sum()); self.nleft = int((~good).sum())
        vp = -b
        with np.errstate(divide='ignore', invalid='ignore'):
            logG = (Delta*np.log(4*fd['rs']) + np.log(np.abs(Sr)) + 2*vp*np.log(fd['E0']/fd['A0']) + (D1 + D2)*np.log(fd['A0']))
        logG[~conv] = -np.inf
        self.logG = np.full(g['S'].shape, -np.inf); self.logG[ok] = logG
        self.sgn = np.ones(g['S'].shape); self.sgn[ok] = np.sign(Sr)

    def moment(self, k, nu):
        g = self.g; U, S = g['U'], g['S']
        F = Fnu_np(nu, [0], -np.log(g['us']))[0].real
        rho, rhob = S**2, (S*U)**2
        with np.errstate(divide='ignore', invalid='ignore'):
            logw = np.log(2.0) + 3*np.log(S) + np.log(U) + np.log(rho - rhob) + 0.5*(k - 4)*np.log(rho*rhob) + self.logG
        logw[~g['ok']] = -np.inf
        sh = np.max(logw)
        val = np.sum(np.exp(logw - sh)*self.sgn*F[:, None]*g['W'])
        return val*mp.exp(mp.mpf(sh))

if __name__ == "__main__":
    quick = len(sys.argv) > 1 and sys.argv[1] == 'quick'
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    pr("# scalar validation of the corner completion (see docstring).  R = exact/predicted.")
    pr("#  R_rad: radial series alone with N levels (as e2e/t3);  R_rad2: the same with 2N;  R_lc: radial (converged points) + LC (rest), N levels;")
    pr("#  R_lc': the same with a second LC setting (N, I, Q, M) = (8000, 1200, 6000, 8) instead of (10000, 1500, 8000, 10).")
    for (D1, D2) in ((3.0, 1.75), (1.75, 0.5)):
        pr(f"# D1 = {D1}, D2 = {D2}")
        pr("#  Delta      J     m     beta  k  nu      R_rad        R_rad2       R_lc          R_lc'        |R_lc-R_lc'|  m^2(R_lc-1)   LC pts / left")
        for beta in ((1.0,) if quick else (0.4, 1.0, 2.0)):
            for (Delta, J) in cases_fixed_beta(beta, [10, 20, 40] if quick else ([10, 20, 40, 70, 100] if beta < 1.5 else [10, 20, 40, 70])):
                t0 = time.time()
                m, b_, eta = info(Delta, J)
                N = int(min(10000, max(4000, 400*Delta)))
                mr = ScalarLC(Delta, J, D1, D2, N, use_lc=False)
                mr2 = ScalarLC(Delta, J, D1, D2, min(2*N, 20000), use_lc=False)
                ml = ScalarLC(Delta, J, D1, D2, N)
                ml2 = ScalarLC(Delta, J, D1, D2, N, lcpar=dict(N=8000, I=1200, Q=6000, M=8))
                for k in (2, 3):
                    for nu in (0.0, 3.0):
                        pd = prediction(Delta, J, D1, D2, k, nu).real
                        R = [mr.moment(k, nu)/pd, mr2.moment(k, nu)/pd, ml.moment(k, nu)/pd, ml2.moment(k, nu)/pd]
                        pr(f"  {Delta:9.4f} {J:4d} {float(m):7.2f} {float(b_):5.3f} {k} {nu:4.1f}  " + "  ".join(f"{float(x):.9f}" for x in R)
                           + f"  {float(abs(R[2]-R[3])):.1e}   {float(m*m*(R[2]-1)):10.4f}    {ml.nlc} / {ml.nleft}")
                pr(f"#    ({time.time()-t0:.0f} s)")
    if not quick:
        open(os.path.join(HERE, 'out_v1_scalar.txt'), 'w').write("\n".join(out) + "\n")
