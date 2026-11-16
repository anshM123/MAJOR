"""Referee: independent reproduction of -C_{2,nu}|_grav = 8 pi G/(nu^2+9/4), d = 3 (numpy/scipy; no code shared with adslight/).
Independent ingredients:
 * T-pole residue SHAPE from the Costa-Goncalves-Penedones leading-Regge Mack coefficient
   ~ 1/(m! (Delta-h+1)_m Gamma(Df-(Delta-J)/2-m)^2), Delta=3, J=2, h=3/2  ->  Gamma(a+m)^2/(m! Gamma(5/2+m)), a = 3/2-Df;
 * overall normalisation NOT CHMRS's R_m constant: fixed by the flat limit, sum_m R_m = -8 pi G c_1 (Penedones, degree-1 amplitude
   8 pi G s^2/(-u)), contacts normalised so that s^2+t^2+u^2 <-> -C_2 = 2 (Barnes lemma); c_1/c_2 = 1/(4 Df - 1);
 * only the kernel Gamma(z)^2 Gamma(c-z) Gamma(cbar-z), c = Df + 1/4 - i nu/2, is CHMRS's.
 * contour Re z = Df - 1/4, trapezoid rule in Im z (exponentially convergent), pole sum to m = 2e5 + power-law tail extrapolation.
"""
import numpy as np
from scipy.special import loggamma, gammaln

def minusC2(Df, nu, G=1/(8*np.pi), M=200000, h=0.01, Y=35.0):
    a = 1.5 - Df
    c = Df + 0.25 - 0.5j*nu
    m = np.arange(M)
    shape = np.exp(2*gammaln(a+m) - gammaln(m+1) - gammaln(2.5+m))   # Gamma(a+m)^2 > 0
    S0 = np.exp(2*gammaln(a) + gammaln(2.5-2*a) - 2*gammaln(2.5-a))   # Gauss 2F1(a,a;5/2;1)
    y = np.arange(-Y, Y+h/2, h)
    z = Df - 0.25 + 1j*y
    U = 2*Df - 2*z
    At = np.zeros_like(z)
    for i0 in range(0, M, 5000):
        mm = m[i0:i0+5000]; sh = shape[i0:i0+5000]
        At += (sh[None, :]/(U[:, None] - 1 - 2*mm[None, :])).sum(1)
    # tail m >= M: shape ~ C m^{2a-7/2}, 1/(U-1-2m) ~ -1/(2m)  ->  -C/2 * M^{2a-9/2+1}/(-(2a-9/2+1))
    p = 2*a - 3.5
    Cm = shape[-1]/(M-1)**p
    At += -0.5*Cm*M**p/(-p)   # sum_{m>=M} C m^p (-1/(2m)) ~ -C/2 * M^p/(-p)  (~1e-15, negligible)
    kern = np.exp(2*loggamma(z) + loggamma(c - z) + loggamma(np.conj(c) - z))
    I = (kern*At).sum()*h/(2*np.pi)
    pref = np.exp(gammaln(2*Df+0.5) - 4*np.real(loggamma(c)))
    return pref*(-8*np.pi*G/(4*Df-1)/S0)*I

for Df in [1.2, 1.75, 2.6]:
    for nu in [0.0, 1.0, 3.0, 7.0, 12.0]:
        v = minusC2(Df, nu)
        t = 1/(nu**2 + 2.25)
        print('Df=%.2f nu=%5.1f  -C2/(8piG) = %.13f  1/(nu^2+9/4) = %.13f  rel=%.1e   imag=%.1e' % (Df, nu, v.real, t, abs(v.real/t-1), abs(v.imag)), flush=True)
