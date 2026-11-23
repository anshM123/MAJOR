"""
pwimg/pw_common.py (PW_IMAGES.md): shared helpers.  Read-only imports of eestep/ee_blocks.py and eestep/ee_quad.py (sys.path; no
file of theirs is written; sys.dont_write_bytecode).

Objects (EESTEP_CLOSE §1; H7 frame x1 = -r yh, x3 = r yh, x2 = -xh, x4 = xh, |x| = 1, eta = cosh B, rho = r e^B, rhob = r e^-B):
  Gh^{Dp}_b(r, eta; mu) = sum_n r^n sum_m B^{Dp}_b[n, m](mu) P_{l+m}(eta),   l = -1/2 + i mu   (continuous-spin (0, b) block / (4r)^Dp)
  E^{Dp}_mu(x; pols)   = sum_{b=0}^{4} C_b(x; pols) v^{b/2} Gh^{Dp}_b(x; mu),  C_0 = c_0, C_b = c_b + c_-b   (ee_quad.laurent_cb)
  raw EE per-nu kernel: K^EE_nu = N(nu) r^2 E^{1}_nu;   f = v' - u' = eta (1 + r^2)/(2 r)  (nu-independent factor of the raw integrand).
"""
import os, sys
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
sys.path.insert(0, os.path.join(HERE, '..', 'eestep'))
import numpy as np
import mpmath as mp
import ee_blocks as EB      # read-only
import ee_quad as Q         # read-only


def Cb_of(r, beta, pol):
    cb, tail, odd, q = Q.laurent_cb(r, beta, pol)
    return np.array([cb[0]] + [cb[b] + cb[-b] for b in range(1, 5)])


def legendre_row(mu, beta, N, dps=40):
    """P_{l+m}(cosh beta), m = -N..N, complex l = -1/2 + i mu, by mpmath (direct; adequate for N <= 80)."""
    with mp.workdps(dps):
        z = mp.cosh(mp.mpf(beta))
        l = mp.mpf(-0.5) + 1j*mp.mpmathify(mu)
        return np.array([complex(mp.legenp(l + m, 0, z, type=3)) for m in range(-N, N + 1)])


def Gh_val(mu, b, r, beta, N=60, Dp=0.0, Prow=None):
    B, m = EB.block_coeffs_b(mu, b, N, Dp)
    P = legendre_row(mu, beta, N) if Prow is None else Prow
    return np.sum((r**np.arange(N + 1))[:, None]*B*P[None, :])


def E_val(mu, r, beta, pol, N=60, Dp=0.0):
    C = Cb_of(r, beta, pol)
    v = EB.v_cross(r*np.exp(beta), beta)
    P = legendre_row(mu, beta, N)
    return sum(C[b]*v**(b/2)*Gh_val(mu, b, r, beta, N, Dp, P) for b in range(5)), C


def cfac(mu):
    """EE/scalar leading normalisation c(mu) = 32(4mu^2+1)(4mu^2+9)/((4mu^2+25)(4mu^2+81))  (EE_FUNCTIONAL §2)."""
    m2 = 4*mu*mu
    return 32*(m2 + 1)*(m2 + 9)/((m2 + 25)*(m2 + 81))
