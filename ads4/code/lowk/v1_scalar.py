"""
Validation on scalar partial waves (no spinning operators), <phi_1 phi_2 phi_1 phi_2>, s-channel (12) -> (34), k = 2, 3.
  exact:   (-1)^J Pihat_{k,nu}[u'^{-(D1+D2)/2} g^{a,b}]  on a (u, t) Gauss grid, with
           R_rad: the radial series alone on every point (N levels, no convergence mask; as e2e/t3_scalar_mixed.py),
           R_lc:  radial series on its converged points + the LC/tail representation (setting A of corner.py) on the others;
                  error budget as for the spinning case (3|A - B| per LC point; power-law extrapolation of the points where the
                  lightcone sum has not converged), 3 x max(|A-B|, |A-C|) over the three tail settings, on the e2e grid and on the corner-refined grid (corner.UP_REF);
  predicted: the derived scalar heavy formula (SPINNING_E2E.md sect. 4.1; e2e/t3_scalar_mixed.prediction), nothing fitted.
Cases:  (D1, D2) = (3, 1.75): corner weight u^{2 D2 + k - 3/2} = u^4 (k = 2): the radial series alone is converged (consistency);
        (D1, D2) = (3.75, 0.5): corner weight u^{3/2} (k = 2), u^{5/2} (k = 3), the corner weight of the spinning k = 2, 3 moments.
Usage: python3 v1_scalar.py [quick]       Output: out_v1_scalar.txt
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np, mpmath as mp
from moments import make_grid, frame_data, series_converged
from blocks import CBlock
from lcblock import LCBlock
from corner import PAR_A, PAR_B, PAR_C, SAFETY, UP_E2E, UP_REF, TensorMomentsLC
from t3_scalar_mixed import prediction
from cases import cases_fixed_beta, info
mp.mp.dps = 25

class ScalarLC(TensorMomentsLC):
    """duck-typed TensorMomentsLC (one 'mode' (0,0)) so that moment_err / _density can be re-used."""
    def __init__(self, Delta, J, D1, D2, N, use_lc=True, upanels=UP_E2E, itail_max=1e-13):
        self.Delta, self.J, self.l = Delta, J, 0
        self.Dphi = D1 + D2 - 1                    # so that the t-law exponent a = 2 Dphi + 2k + 2l of moment_err is 2 D1 + 2 D2 + 2k - 2
        g = make_grid(Delta, J, upanels=upanels); self.g = g
        ok = g['ok']
        fd = frame_data(g['U'][ok], g['S'][ok])
        b = (D1 - D2)/2
        assert b > 0
        conv = series_converged(fd['rs'], fd['chi'], J, N) if use_lc else np.ones(fd['rs'].size, bool)
        Sr = np.zeros(fd['rs'].size); dS = np.zeros(fd['rs'].size)
        blk = CBlock(Delta, J, b, b, N)
        Sr[conv] = blk.series(fd['rs'][conv], fd['chi'][conv], sig=+1)[:, 0, 0]
        del blk
        lc = np.zeros(fd['rs'].size, bool)
        self.nlc = self.nleft = 0
        if use_lc and (~conv).any():
            idx = np.where(~conv)[0]
            lbA = LCBlock(Delta, J, b, b, lmax=0, **PAR_A)
            SA = lbA.series(fd['rs'][idx], fd['chi'][idx], kmax=0, lmax=0)[:, 0, 0]; tA = lbA.i_tail.copy()
            lbB = LCBlock(Delta, J, b, b, lmax=0, coeffs=lbA.C, **PAR_B)
            SB = lbB.series(fd['rs'][idx], fd['chi'][idx], kmax=0, lmax=0)[:, 0, 0]; tB = lbB.i_tail.copy()
            lbC = LCBlock(Delta, J, b, b, lmax=0, coeffs=lbA.C, **PAR_C)
            SC = lbC.series(fd['rs'][idx], fd['chi'][idx], kmax=0, lmax=0)[:, 0, 0]; tC = lbC.i_tail.copy()
            good = (tA < itail_max) & (tB < itail_max) & (tC < itail_max)
            Sr[idx[good]] = SA[good]; dS[idx[good]] = np.maximum(np.abs(SA - SB), np.abs(SA - SC))[good]/np.abs(SA[good])
            conv[idx[good]] = True; lc[idx[good]] = True
            self.nlc = int(good.sum()); self.nleft = int((~good).sum())
        vp = -b
        with np.errstate(divide='ignore', invalid='ignore'):
            logG = (Delta*np.log(4*fd['rs']) + np.log(np.abs(Sr)) + 2*vp*np.log(fd['E0']/fd['A0']) + (D1 + D2)*np.log(fd['A0']))
        logG[~conv] = -np.inf
        sh = g['S'].shape
        self.logscale = np.full(sh, -np.inf); self.logscale[ok] = logG
        md = np.zeros(sh); md[ok] = np.where(conv, np.sign(Sr), 0.0)
        dm = np.zeros(sh); dm[ok] = dS
        self.modes = {(0, 0): md}; self.dmodes = {(0, 0): dm}
        self.conv = np.zeros(sh, bool); self.conv[ok] = conv
        self.lcpts = np.zeros(sh, bool); self.lcpts[ok] = lc

    def value(self, k, nu):
        (v, sh), err_lc, err_drop, share = self.moment_err(k, nu, 0, 0)
        return v*mp.exp(mp.mpf(sh)), err_lc, err_drop, share

if __name__ == "__main__":
    quick = len(sys.argv) > 1 and sys.argv[1] == 'quick'
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    pr("# scalar validation of the corner completion (see docstring).  R = exact/predicted;  err = err_lc + err_drop (relative);")
    pr("# R_lc(e2e grid), R_lc(refined grid) and their difference (u-quadrature estimate);  share = part of the moment from LC points.")
    for (D1, D2) in ((3.75, 0.5), (3.0, 1.75)):
        pr(f"# D1 = {D1}, D2 = {D2}")
        pr("#  Delta      J     m     beta  k  nu      R_rad         R_lc(e2e grid)  err      R_lc(refined)   err      |diff|   share    m^2(R_lc,ref - 1)")
        for beta in ((1.0,) if quick else (0.4, 1.0)):
            for (Delta, J) in cases_fixed_beta(beta, [10, 20, 40] if quick else [10, 20, 40, 70, 100]):
                t0 = time.time()
                m, b_, eta = info(Delta, J)
                N = int(min(10000, max(4000, 400*Delta)))
                mr = ScalarLC(Delta, J, D1, D2, N, use_lc=False)
                ml = ScalarLC(Delta, J, D1, D2, N)
                mf = ScalarLC(Delta, J, D1, D2, N, upanels=UP_REF)
                for k in (2, 3):
                    for nu in (0.0, 3.0):
                        pd = prediction(Delta, J, D1, D2, k, nu).real
                        vr = mr.value(k, nu)[0]
                        vl, el, ed, sl = ml.value(k, nu); vf, fl, fd_, sf = mf.value(k, nu)
                        Rr, Rl, Rf = mp.re(vr)/pd, mp.re(vl)/pd, mp.re(vf)/pd
                        pr(f"  {Delta:9.4f} {J:4d} {float(m):7.2f} {float(b_):5.3f} {k} {nu:4.1f}  {float(Rr):.9f}  {float(Rl):.9f}  {el+ed:.0e}  "
                           f"{float(Rf):.9f}  {fl+fd_:.0e}  {float(abs(Rl-Rf)):.0e}  {sf:.0e}   {float(m*m*(Rf-1)):10.4f}")
                pr(f"#    ({time.time()-t0:.0f} s; LC points / left: e2e grid {ml.nlc}/{ml.nleft}, refined {mf.nlc}/{mf.nleft})")
    if not quick:
        open(os.path.join(HERE, 'out_v1_scalar.txt'), 'w').write("\n".join(out) + "\n")
