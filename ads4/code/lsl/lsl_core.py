"""
lsl/lsl_core.py (LARGESPIN_L.md): shared helpers, re-implemented from the published definitions.

  * NF02 weights psi_2, psi_3: ads4/prl/data/nf02_psi.json (exact rationals, polynomials of degree 14, psi(0) = 0, (1-p)^5 | psi).
  * Phi_k(x) = int_0^1 psi_k(p) cos(p x) dp in CLOSED FORM (repeated integration by parts of a polynomial; mpmath for small x):
        int_0^1 psi e^{ipx} dp = sum_j (-1)^j [psi^(j)(1) e^{ix} - psi^(j)(0)] / (ix)^{j+1}.
  * Omega_k(beta) = int_0^1 psi_k(p) P_{-1/2+i Lam p}(cosh beta) dp, by Mehler-Dirichlet in the Lemma-C substitution
        Omega_k = (2/pi) int_0^{pi/2} Phi_k(Lam s)/cosh(s/2) dphi,   s = 2 asinh(sinh(beta/2) sin phi)
    with the closed-form Phi (valid at any beta; ls_core.OmLR's Gauss-Legendre Phi is accurate only for Lam*beta <~ 4000).
  * M1(beta) = int_0^inf 2 nu P_{-1/2+i nu}(cosh beta)/(e^{2 pi nu}+1) dnu   (ROUTEA_LARGE_B §4.1).
  * P_{-1/2+i nu}(cosh beta) and d/dbeta for real nu (same substitution).
No project module is imported here; nothing outside code/lsl/cache is written.
"""
import os, sys, json
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(ADS4, 'prl', 'data')


def lam(Dg):
    return np.sqrt((Dg - 2.0)*(Dg - 1.0))


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def mgap(Dg, J):
    """mass of the twist-edge state Delta = Dg + J (A1': twist >= Dg)."""
    return mphys(Dg + J, J)


def beta_of(J, m):
    return 2*np.arcsinh((J + 0.5)/m)


def bhat_of(J, m, Dg=1000):
    return (2*J + 1)*lam(Dg)/m


_PSI = {}


def psi_coeffs(k):
    """exact rational coefficients c_n of psi_k(p) = sum c_n p^n (k0 included)."""
    if k not in _PSI:
        d = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
        k0 = Fr(d['k0'])
        _PSI[k] = [Fr(c)/k0 for c in d['psi%d_exact' % k]]
    return _PSI[k]


def psi(k, p):
    c = [float(x) for x in psi_coeffs(k)]
    return np.polynomial.polynomial.polyval(np.asarray(p, float), c)


def derivs_at(k, p0):
    """exact psi^(j)(p0), j = 0..deg (p0 = 0 or 1)."""
    c = psi_coeffs(k); D = len(c) - 1
    out = []
    for j in range(D + 1):
        s = Fr(0)
        for n in range(j, D + 1):
            ff = 1
            for i in range(j):
                ff *= (n - i)
            s += c[n]*ff*(Fr(p0)**(n - j))
        out.append(s)
    return out


_D0 = {}; _D1 = {}


def d0(k):
    if k not in _D0:
        _D0[k] = [float(x) for x in derivs_at(k, 0)]
    return _D0[k]


def d1(k):
    if k not in _D1:
        _D1[k] = [float(x) for x in derivs_at(k, 1)]
    return _D1[k]


def Phi(k, x):
    """closed form, vectorised; for x < 40 an mpmath evaluation at 50 digits is used (cancellation)."""
    x = np.atleast_1d(np.asarray(x, float))
    out = np.empty_like(x)
    big = x >= 40
    if big.any():
        xb = x[big]
        a0 = d0(k); a1 = d1(k)
        acc = np.zeros_like(xb, dtype=complex)
        ex = np.exp(1j*xb)
        ix = 1j*xb
        pw = ix.copy()
        for j in range(len(a0)):
            acc += (-1)**j*(a1[j]*ex - a0[j])/pw
            pw = pw*ix
        out[big] = acc.real
    small = ~big
    if small.any():
        out[small] = [_Phi_mp(k, xi) for xi in x[small]]
    return out


def _Phi_mp(k, x):
    with mp.workdps(50):
        c = psi_coeffs(k)
        if x == 0:
            return float(sum(mp.mpf(cn.numerator)/cn.denominator/(n + 1) for n, cn in enumerate(c)))
        xm = mp.mpf(x)
        a0 = [mp.mpf(v.numerator)/v.denominator for v in derivs_at(k, 0)]
        a1 = [mp.mpf(v.numerator)/v.denominator for v in derivs_at(k, 1)]
        ex = mp.expj(xm); ix = 1j*xm; pw = ix; acc = mp.mpc(0)
        for j in range(len(a0)):
            acc += (-1)**j*(a1[j]*ex - a0[j])/pw
            pw *= ix
        return float(acc.real)


def gl01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return (x + 1)/2, w/2


def Omega(k, beta, Dg=1000, nphi=6000):
    """int_0^1 psi_k(p) P_{-1/2+i Lam p}(cosh beta) dp (Mehler-Dirichlet, closed-form Phi)."""
    L = lam(Dg)
    ph, w = gl01(nphi); ph = ph*np.pi/2; w = w*np.pi/2
    s = 2*np.arcsinh(np.sinh(beta/2)*np.sin(ph))
    return (2/np.pi)*(Phi(k, L*s)/np.cosh(s/2)*w).sum()


def P_real(beta, nu, nphi=3000):
    """P_{-1/2+i nu}(cosh beta) and dP/dbeta for an array of real nu."""
    ph, w = gl01(nphi); ph = ph*np.pi/2; w = w*np.pi/2
    S = np.sinh(beta/2)*np.sin(ph)
    s = 2*np.arcsinh(S); ch = np.sqrt(1 + S*S)            # cosh(s/2)
    dsdb = np.cosh(beta/2)*np.sin(ph)/ch
    nu = np.atleast_1d(np.asarray(nu, float))
    C = np.cos(np.outer(nu, s)); Sn = np.sin(np.outer(nu, s))
    P = (2/np.pi)*(C/ch*w).sum(1)
    # d/dbeta [cos(nu s)/cosh(s/2)] = [-nu sin(nu s) - cos(nu s) tanh(s/2)/2] ds/dbeta / cosh(s/2)
    th = S/ch
    dP = (2/np.pi)*(((-nu[:, None]*Sn - C*th/2)*dsdb/ch)*w).sum(1)
    return P, dP


def M1(beta):
    """int_0^inf 2 nu P_{-1/2+i nu}(cosh beta)/(e^{2 pi nu}+1) dnu."""
    nu, w = gl01(400); nu = nu*12; w = w*12
    P, _ = P_real(beta, nu, nphi=max(600, int(60*beta) + 600))
    return (2*nu*P/(np.exp(2*np.pi*nu) + 1)*w).sum()


SIGMA = {'++': 1.0, '+-': -1.0}


def LR_action(beta, u, sig, Dg=1000):
    """leading-Regge per-state action 2 Omega_2 + u sigma Omega_3 (= the hybrid's flat-Regge + curvature)."""
    return 2*Omega(2, beta, Dg) + u*sig*Omega(3, beta, Dg)
