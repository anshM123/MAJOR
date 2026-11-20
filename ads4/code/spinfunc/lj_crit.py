"""
lj_crit (HLJ_ANALYTIC.md §3.4): search for OTHER critical points of the travelling-wave exponent Phi_+- (lj_saddle.phi) near the
scattering crystal, at given (Delta, J, nu).  A non-oscillating contribution to H_+- would need a critical point z' with
|Re kappa~(z')| << beta (local nu-frequency ~ 0) and Re Phi(z') not far below Re Phi(z_+-).
Random starts: z = z_c + a * (sigma_b * N(0,1) complex in the 6 base coordinates) and x~0 ~ U(-0.95, 0.95) * (1 + 0.3 i N(0,1)).
Usage: python3 lj_crit.py Delta J nu [nstart] [sigma_b]  -> prints the distinct critical points, sorted by Re Phi.
"""
import os, sys
os.environ.setdefault('L1_DEG', '2')
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l1_tps as T
import l1_crystal as Cr
import lj_saddle as S


def newton_free(z, nu, Delta, J, sgn, it=80):
    for k in range(it):
        zv = [T.var(i, z[i]) for i in range(T.NV)]
        try:
            P, _, _, _ = S.phi(zv, nu, Delta, J, sgn)
            g, H = T.grad_hess(P)
            dz = np.linalg.solve(H, g)
        except Exception:
            return z, False
        if not np.all(np.isfinite(dz)):
            return z, False
        sc = max(1e-3, np.abs(z[1:5]).max())
        lam = min(1.0, 0.3*sc/max(np.abs(dz[1:]).max(), 1e-300), 0.3/max(abs(dz[0]), 1e-300))
        z = z - lam*dz
        if np.abs(dz).max() < 1e-11*max(1.0, np.abs(z).max()):
            return z, True
    return z, False


def search(Delta, J, nu, sgn=1, nstart=60, sig=0.5, seed=1):
    rng = np.random.default_rng(seed)
    zc, a = Cr.saddle(nu, Delta)
    zp, ok, _ = S.track_J(nu, Delta, J, sgn)
    Ip = S.info(zp, nu, Delta, J, sgn)
    found = [(zp, Ip)]
    for n in range(nstart):
        z = zc.copy()
        z[1:] = zc[1:] + a*sig*(rng.normal(size=6) + 1j*rng.normal(size=6))
        z[5:] = a*sig*(rng.normal(size=2) + 1j*rng.normal(size=2))
        z[0] = rng.uniform(-0.95, 0.95)*(1 + 0.3j*rng.normal())
        z, ok = newton_free(z, nu, Delta, J, sgn)
        if not ok:
            continue
        I = S.info(z, nu, Delta, J, sgn)
        if I['gnorm'] > 1e-6 or not np.isfinite(I['Phi']):
            continue
        if any(np.abs(z - f[0]).max() < 1e-7*max(1, np.abs(z).max()) for f in found):
            continue
        found.append((z, I))
    return Ip, found


if __name__ == '__main__':
    Delta = float(sys.argv[1]); J = int(sys.argv[2]); nu = float(sys.argv[3])
    nstart = int(sys.argv[4]) if len(sys.argv) > 4 else 60
    sig = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
    alpha = Delta - J - 2.0; m = np.sqrt(alpha*(alpha + 2*J + 1)); beta = np.arccosh(1 + (2*J + 1)**2/(2*m*m))
    for sgn in (1, -1):
        Ip, found = search(Delta, J, nu, sgn, nstart, sig)
        print('Delta=%g J=%d nu=%g wave %+d: %d distinct critical points (beta=%.4f)' % (Delta, J, nu, sgn, len(found), beta))
        for z, I in sorted(found, key=lambda f: -f[1]['Phi'].real):
            kt = I['kappa'] - np.log(2)
            print('   Re(Phi-Phi_s)=%+10.3f  Re kappa~/beta=%+8.4f  Im kappa~=%+.4f  t=%s  x0=%s  min|eig|=%.3g' % (
                (I['Phi'] - Ip['Phi']).real, kt.real/beta, kt.imag, np.round(I['t'], 6), np.round(z[0], 4), I['mineig']))
