r"""
lrc3_ee (LRC_CLOSE.md §4): what carries over from the scalar j'=0 kernel (lrc1, lrc2) to the EE kernels of C'_EE and C4^ds.

(a) The EE (2,4) dressing.  Both C'_EE (RESULT_D_STATEMENT §2) and C4^ds (RESULT_Q_STAGE12 §1.1) use
        EE = [V_{0,42}^2 (H24 - 2 V_{2,04} V_{4,02})/(X20 X40)]^2   together with   V_{0,42}^{-4}   (repro_engine/integrand.py),
    so the EE kernel is the scalar kernel's (theta, z0) integrand times  D = (H24 - 2 V_{2,04} V_{4,02})^2/(X20 X40)^2 (nu- and z0-free).
    D contains no z0: checked on the (theta, phi) grid for random null z2, z4.  Its x0 (theta) dependence is NOT small: it is a bounded
    nu-independent factor, of relative size up to ~ (1 - rho)^-8 (denominators X20, X40, X04 ~ -(1 -+ rho cos theta)): printed.
(b) Pole content of the leading light-cone (collinear) coefficients  (h-a)_p (h+b)_p/((2h)_p p!),  h = (Dp + l)/2 = (1/2 + i nu)/2:
    with a = 0 (the 1<->3-symmetric light-ray pair) the poles in the de-shift half-strip 0 < Im nu <= 3 sit only at Im nu = 3/2, for every b;
    with a = 1/2 (control: a 1<->3-asymmetric pair) a pole appears at Im nu = 1/2, inside the k = 2 strip.
Output: cache/lrc3_ee.txt.  seconds.
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
import mpmath as mp

OUT = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def dot(u, v):
    return -u[..., 0]*v[..., 0] + u[..., 1]*v[..., 1] + u[..., 2]*v[..., 2]


def Vf(xi, z, xj, xk):
    """V_{i,jk}(z) = (z.(xj-xi) (xi-xk)^2 - z.(xk-xi) (xi-xj)^2)/(xj-xk)^2  (as repro_engine/integrand.py)."""
    return (dot(z, xj - xi)*dot(xi - xk, xi - xk) - dot(z, xk - xi)*dot(xi - xj, xi - xj))/dot(xj - xk, xj - xk)


def Hf(xi, zi, xj, zj):
    return dot(zi, zj)*dot(xi - xj, xi - xj) + 2*dot(zi, xj - xi)*dot(zj, xi - xj)


def null_vec(rng):
    a = rng.uniform(0, 2*np.pi)
    return np.array([1.0, np.cos(a), np.sin(a)])*rng.uniform(0.5, 2.0)


def part_a():
    pr('(a) EE dressing D = EE * V_{0,42}^{-4} = (H24 - 2 V_{2,04} V_{4,02})^2/(X20 X40)^2: z0-independent, nu-independent, bounded (not small) in x0')
    rng = np.random.default_rng(7)
    th = 2*np.pi*np.arange(96)/96; ph = 2*np.pi*np.arange(97)/97
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    for (beta, rho) in ((1.0, 0.3), (5.0, 0.3), (5.0, 0.8), (12.0, 0.8), (12.0, 0.5)):
        r = rho*np.exp(-beta)
        # everything in the rest frame of xhat (D is Lorentz invariant; this avoids e^beta cancellations): x2 = -yhat, x4 = yhat,
        # x0 = r Lambda_beta^{-1} n_theta (components ~ rho), z0 = Lambda_beta^{-1}(1, cos phi, sin phi), polarisations O(1) null vectors.
        Li = np.array([[np.cosh(beta), -np.sinh(beta), 0], [-np.sinh(beta), np.cosh(beta), 0], [0, 0, 1.0]])
        yh = np.array([1.0, 0.0, 0.0])
        x2 = np.broadcast_to(-yh, TH.shape + (3,)); x4 = np.broadcast_to(yh, TH.shape + (3,))
        x0 = r*np.einsum('ij,abj->abi', Li, np.stack([0*TH, np.cos(TH), np.sin(TH)], -1))
        z0 = np.einsum('ij,abj->abi', Li, np.stack([np.ones_like(PH), np.cos(PH), np.sin(PH)], -1))
        wz = []; wx = []
        for _ in range(4):
            z2 = np.broadcast_to(null_vec(rng), TH.shape + (3,)); z4 = np.broadcast_to(null_vec(rng), TH.shape + (3,))
            V042 = Vf(x0, z0, x4, x2)
            E = V042**2*(Hf(x2, z2, x4, z4) - 2*Vf(x2, z2, x0, x4)*Vf(x4, z4, x0, x2))/(dot(x2 - x0, x2 - x0)*dot(x4 - x0, x4 - x0))
            D = E**2/V042**4
            wz.append(np.max(np.abs(D - D[:, :1])/np.abs(D[:, :1])))
            x00 = 0*x0
            D0 = ((Hf(x2, z2, x4, z4) - 2*Vf(x2, z2, x00, x4)*Vf(x4, z4, x00, x2))/(dot(x2, x2)*dot(x4, x4)))**2
            wx.append(np.max(np.abs(D - D0))/np.max(np.abs(D0)))
        pr(f'    beta = {beta:4.1f}, rho = {rho}, r = {r:.1e}: max_z0 rel variation of D = {max(wz):.1e};  '
           f'max|D - D_LR|/max|D_LR| over x0, 4 random polarisation pairs = ' + ', '.join(f'{x:.3f}' for x in wx))


def part_b():
    pr('(b) poles of the collinear coefficients (h-a)_p(h+b)_p/((2h)_p p!), h = (1/2 + i nu)/2, on nu = i y, 0 < y <= 3 (p <= 12)')

    def poles(a, b):
        found = []
        for k in range(1, 61):
            y = mp.mpf(k)/20
            h = (mp.mpf(0.5) + 1j*(1j*y))/2
            for dy in (mp.mpf('1e-9'),):
                hh = (mp.mpf(0.5) + 1j*(1j*(y + dy)))/2
                c = max(abs(mp.rf(hh - a, p)*mp.rf(hh + b, p)/(mp.rf(2*hh, p)*mp.factorial(p))) for p in range(1, 13))
                if c > 1e6:
                    found.append(float(y))
        return found
    for (a, b) in ((0, 0), (0, 0.5), (0, 1), (0, -1), (0, 2), (0, -2), (1, 0), (0.5, 0.5), (-1, 1)):
        pr(f'    a = {a:+}, b = {b:+}: poles at Im nu = {poles(mp.mpf(a), mp.mpf(b))}')


if __name__ == '__main__':
    pr('# lrc3_ee: carrying LRC* from the scalar j\'=0 kernel to the EE kernels (LRC_CLOSE.md §4)')
    part_a(); part_b()
    with open(os.path.join(HERE, 'cache', 'lrc3_ee.txt'), 'w') as f:
        f.write('\n'.join(OUT) + '\n')
