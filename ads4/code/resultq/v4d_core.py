"""
v4d_core (RESULT_Q_V4DPRIME.md): closed forms and exact (Fourier) evaluation of the corner pairings of the CLS-D spinning k = 4 kernel in
the complete gauge of v4c_core (x0 = 0, x2 = -e, x4 = +e, z0 fixed), plus the heavy normalisation and the stored de-shift weights.

(1) u-channel identity W_u = <T4T2><T3T1>, edge residue of t1^3 t3^3 f4 W_u (EE dual structure, T42 = -1/64), timelike tube
    x1 -> x1 + i sgn eps e.  EXACT closed form (derived here; checked against v4c_wu.U_eps to the FFT precision 1e-7):
        X13 = -(1-c)/2 + 2 i sgn eps + eps^2,   N = (z1.z3) X13 - 2 (z1.x13)(z3.x13) = eps^2 (1 + c)    (t-independent),
        Res f4 = -(1 - X13^2)/(4 X13),   U(dphi) = (pi^2/64) T42 eps^4 (1+c)^2 (1 - X13^2)/X13^6,    c = cos(dphi).
(2) Celestial factor:  int dphi1 dphi3 C(phi1, phi3) g(phi1 - phi3) = (2 pi)^2 sum_n fh_n^2 gh_n,  fh_n the Fourier coefficients of
    (2 - 2 cos phi)^s, s = -delta/2:  fh_n = (-1)^n Gamma(2s+1)/(Gamma(s+n+1) Gamma(s-n+1))  (exact),  gh_n = FFT of
    (2 - 2cos dphi)^{(2+delta)/2} g(dphi).  Valid for any nu (no oscillatory quadrature).
(3) Small-dphi expansion of K_nu(dphi) = int dphibar (A10 A30)^{-delta/2} = a(nu) + b(nu) dphi^{-2 i nu} + O(dphi^2):
        a = 2 pi Gamma(1 - 2 delta)/Gamma(1 - delta)^2,   b = 2^{1 + 2 i nu} [B(1/2, 1-delta) + B(delta - 1/2, 1-delta)]/2.
    Leading small-eps pairing (dphi = sqrt(eps) s):
        G_nu(eps) = (pi^2/16) T42 eps^{-1/4} [a J(5/2 + i nu) eps^{i nu/2} + b J(5/2 - i nu) eps^{-i nu/2}] (1 + O(eps)),
        J(sig) = 2^{sig+1} c^{al-6} Gamma(al) Gamma(6-al)/Gamma(6),  al = (sig+1)/2,  c = 2 e^{-i sgn pi/2}.
(4) Heavy normalisation (Stage 4: h_j = N_com R_j, N_com absorbed in chi):  N_com(nu) = 2 pi K_A^ref K_B^ref = 2^11 pi^2 |Gamma(11/4 + i nu/2)|^4
    (K(c, lam) of v3_regge_kernel with c_A = 4 - delta, lam_A = 1 - delta; c_B = 3 + delta, lam_B = delta).  The functional is
    int chi(nu) Psi_nu/N_com(nu) dnu, so the weight on the v4c-normalised position kernel is w(nu) = chi(nu) N_com(0)/N_com(nu) (unit:
    the v4c kernel at nu = 0).
"""
import os, sys, json
from fractions import Fraction
import numpy as np
import mpmath as mp
from scipy.special import loggamma

HERE = os.path.dirname(os.path.abspath(__file__))
T42 = -1.0/64


def U_exact(dph, eps, sgn=+1):
    c = np.cos(dph)
    X = -(1 - c)/2 + 2j*sgn*eps + eps**2
    return (np.pi**2/64)*T42*eps**4*(1 + c)**2*(1 - X**2)/X**6


def fhat(nu, nmax):
    """Fourier coefficients fh_n, n = 0..nmax, of (2 - 2 cos phi)^s, s = -(1/2 + i nu)/2 (complex128 recurrence; fh_0 by mpmath)."""
    s = -(0.5 + 1j*nu)/2
    f0 = complex(mp.gamma(2*mp.mpc(s) + 1)/mp.gamma(mp.mpc(s) + 1)**2)
    n = np.arange(nmax)
    r = -(s - n)/(s + n + 1)
    return f0*np.concatenate([[1.0 + 0j], np.cumprod(r)])


def pair_exact(nu, gfun, N=2**16):
    """(2 pi)^2 sum_n fh_n^2 gh_n with gh_n the Fourier coefficients of A13^{(2+delta)/2} gfun(dphi) on [-pi, pi)."""
    d = 0.5 + 1j*nu
    x = -np.pi + 2*np.pi*np.arange(N)/N
    A13 = 2 - 2*np.cos(x)
    with np.errstate(divide='ignore', invalid='ignore'):
        g = np.where(A13 > 0, A13**((2 + d)/2), 0)*gfun(x)
    gh = np.fft.fft(g)/N*np.exp(1j*np.pi*np.fft.fftfreq(N, 1.0/N))      # start point -pi: e^{-i n (-pi)}
    nn = np.fft.fftfreq(N, 1.0/N).astype(int)
    fh = fhat(nu, N//2 + 1)
    F = fh[np.abs(nn)]
    return (2*np.pi)**2*np.sum(F**2*gh)


def ab(nu):
    d = mp.mpf(1)/2 + 1j*mp.mpf(nu)
    a = 2*mp.pi*mp.gamma(1 - 2*d)/mp.gamma(1 - d)**2
    I = (mp.beta(mp.mpf(1)/2, 1 - d) + mp.beta(d - mp.mpf(1)/2, 1 - d))/2
    return a, 2**(1 + 2j*mp.mpf(nu))*I


def Jfun(sig, sgn=+1):
    al = (sig + 1)/2
    c = 2*mp.exp(-sgn*1j*mp.pi/2)
    return 2**(sig + 1)*c**(al - 6)*mp.gamma(al)*mp.gamma(6 - al)/mp.gamma(6)


def coef_lead(nu, sgn=+1):
    """(cp, cm): G_nu(eps) ~ (pi^2/16) T42 eps^{-1/4} [cp eps^{i nu/2} + cm eps^{-i nu/2}]."""
    a, b = ab(nu)
    nu = mp.mpf(nu)
    return complex(a*Jfun(mp.mpf(5)/2 + 1j*nu, sgn)), complex(b*Jfun(mp.mpf(5)/2 - 1j*nu, sgn))


def G_lead(nu, eps, sgn=+1, cpm=None):
    cp, cm = cpm if cpm is not None else coef_lead(nu, sgn)
    return np.pi**2/16*T42*eps**(-0.25)*(cp*eps**(0.5j*nu) + cm*eps**(-0.5j*nu))


def logNcom(nu):
    nu = np.asarray(nu, float)
    return np.log(2.0**11*np.pi**2) + 4*loggamma(11/4 + 0.5j*nu).real


def weights():
    """{name: (nu, chi, dnu)} for the stored (Stage 4) and extended-precision (Stage 7) C4 weights (value, derivative)."""
    out = {}
    s = json.load(open(os.path.join(HERE, 'functionals_q4_deshift.json')))['EE3']
    for nm, key in (('stored_value', 'chi_value'), ('stored_deriv', 'chi_deriv')):
        chi = np.array([float(v) for v in s[key]]); dnu = float(s['dnu'])
        out[nm] = ((np.arange(len(chi)) + 0.5)*dnu, chi, dnu)
    x = json.load(open(os.path.join(HERE, 'functionals_q7_xprec.json')))
    for nm, key in (('xp_value', 'C4_value'), ('xp_deriv', 'C4_deriv')):
        chi = np.array([float(v) for v in x[key]['chi']]); dnu = float(Fraction(x[key]['dnu']))
        out[nm] = ((np.arange(len(chi)) + 0.5)*dnu, chi, dnu)
    return out


def weights_mp(name, dps=70):
    """XP weights at full (64-digit) precision as mpf list, with nu grid (mpf)."""
    mp.mp.dps = dps
    x = json.load(open(os.path.join(HERE, 'functionals_q7_xprec.json')))
    key = {'xp_value': 'C4_value', 'xp_deriv': 'C4_deriv'}[name]
    dnu = mp.mpf(Fraction(x[key]['dnu']).numerator)/Fraction(x[key]['dnu']).denominator
    chi = [mp.mpf(v) for v in x[key]['chi']]
    nu = [(k + mp.mpf(1)/2)*dnu for k in range(len(chi))]
    return nu, chi, dnu


def RtEE(c):
    """closed-form [TT]_0 edge residue (EE, k = 4) of v4c_mft (includes (-2 pi i)^2)."""
    return -4*np.pi**2*(c - 3)*(c - 1)*(c + 1)/524288.0
