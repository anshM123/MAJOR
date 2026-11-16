"""
Formulas used for the predictions (d = 3 throughout).  Functions copied/adapted from reggemoment/chmrs_formulas.py
(CHMRS = arXiv:2106.10274):
  (2.28)-(2.30)  alpha = Delta-J-2, alphabar = Delta+J-1, m^2 = alpha alphabar, cosh beta = eta_AdS = (alpha^2+alphabar^2)/(2 m^2)
  (2.55)         gamma_a(nu) = |Gamma((a - 1/2 + i nu)/2)|^2
  (3.6)          b_{Delta,J}; generalised here to Delta_34 != 0 (b_coeff)
  harmonic functions  F^nu_mu(xi) = (1/2pi) int_0^{2pi} dtheta (cosh xi - sinh xi cos theta)^{-1/2+i nu} cos(mu theta)
                      F^nu_0 = P_{-1/2+i nu}(cosh xi).
"""
import mpmath as mp, numpy as np

def alpha_pair(Delta, J):
    return Delta - J - 2, Delta + J - 1

def m2_eta(Delta, J):
    a, ab = alpha_pair(Delta, J)
    return a*ab, (a*a + ab*ab)/(2*a*ab)

def DeltaJ_from(m, beta_J):
    """inverse map for integer J: returns Delta with given J and m (alpha = m e^{-beta/2}...) is not used; see cases()."""
    raise NotImplementedError

def gamma_a(a, nu):
    return mp.gamma((a - mp.mpf(1)/2 - 1j*nu)/2)*mp.gamma((a - mp.mpf(1)/2 + 1j*nu)/2)

def log_b(Delta, J, D34=0):
    """log of the Lorentzian shadow coefficient, d = 3;  D34 = 0: CHMRS (3.6).  For D34 != 0 the squares
    Gamma((J+Delta)/2)^2, Gamma((Delta-1-J)/2)^2 are replaced by Gamma(. + D34/2) Gamma(. - D34/2)."""
    g = mp.loggamma
    h = mp.mpf(D34)/2
    return (mp.log(mp.pi) + g(Delta-1) + g(J+mp.mpf(1)/2) + g((J+Delta)/2 + h) + g((J+Delta)/2 - h)
            + g((Delta-1-J)/2 + h) + g((Delta-1-J)/2 - h)
            - mp.log(2) - g(Delta-mp.mpf(1)/2) - g(J+1) - g(J+Delta) - g(Delta-1-J))

def Pnu(nu, eta):
    return mp.hyp2f1(mp.mpf(1)/2 - 1j*nu, mp.mpf(1)/2 + 1j*nu, 1, (1-mp.mpf(eta))/2)

def Fnu_mp(nu, mu, xi):
    """F^nu_mu(xi) by direct quadrature (mpmath)."""
    xi = mp.mpf(xi)
    f = lambda th: (mp.cosh(xi) - mp.sinh(xi)*mp.cos(th))**(-mp.mpf(1)/2 + 1j*nu)*mp.cos(mu*th)
    return mp.quad(f, mp.linspace(0, mp.pi, 9))/mp.pi

def Fnu_closed(nu, mu, xi):
    """closed form: F^nu_mu = Gamma(1/2+i nu)/Gamma(1/2+i nu+|mu|) P^{|mu|}_{-1/2+i nu}(cosh xi)  (Legendre function of type 3)."""
    mu = abs(mu)
    return mp.gamma(mp.mpf(1)/2 + 1j*nu)/mp.gamma(mp.mpf(1)/2 + 1j*nu + mu)*mp.legenp(-mp.mpf(1)/2 + 1j*nu, mu, mp.cosh(mp.mpf(xi)), type=3)

def Fnu_np(nu, mus, xi, M=None):
    """numpy version by the periodic trapezoid rule; xi array; returns dict mu -> complex array."""
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    exi = np.exp(np.max(xi))
    if M is None:
        M = int(2**np.ceil(np.log2(600 + 60*exi + 40*abs(nu)*(1 + np.max(xi)))))
    th = 2*np.pi*(np.arange(M) + 0.5)/M
    out = {mu: np.zeros(xi.size, dtype=complex) for mu in mus}
    ch, sh = np.cosh(xi), np.sinh(xi)
    for i in range(xi.size):
        base = np.exp((-0.5 + 1j*nu)*np.log(ch[i] - sh[i]*np.cos(th)))
        for mu in mus:
            out[mu][i] = np.mean(base*np.cos(mu*th))
    return out
