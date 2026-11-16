"""
Formulas quoted from Caron-Huot, Mazac, Rastelli, Simmons-Duffin, arXiv:2106.10274 ("CHMRS"), transcribed from the
LaTeX source (not from the PDF text).  Everything is for general d; the scripts use d = 3.

  (2.28)-(2.30)  alpha = Delta-J-d+1, alphabar = Delta+J-1, m^2 = alpha*alphabar, eta_AdS = (alpha^2+alphabar^2)/(2 m^2)
  (2.14)         n_J      (D = d+1)
  (2.33)         p^MFT_{Delta,J}
  (2.32)         q_{Delta,J}
  (2.55)         gamma_a(nu)
  (2.57)         a_{Delta_phi}
  (3.6)          b_{Delta,J}
  (B.1)          P_J(eta) = 2F1(-J, J+d-2; (d-1)/2; (1-eta)/2)
  (3.27)         first 1/m^2 correction of C_{k,nu}[Delta,J]
"""
import mpmath as mp

def alpha_pair(Delta, J, d=3):
    return Delta - J - d + 1, Delta + J - 1

def m2_eta(Delta, J, d=3):
    a, ab = alpha_pair(Delta, J, d)
    m2 = a*ab
    return m2, (a*a + ab*ab)/(2*m2)

def n_J(J, d=3):
    D = d + 1
    return 2**D*mp.pi**(mp.mpf(D-2)/2)/mp.gamma(mp.mpf(D-2)/2)*mp.rf(J+1, D-4)*(2*J+D-3)

def log_pMFT(Delta, J, Dphi, d=3):
    g = mp.loggamma
    h = mp.mpf(d)/2
    r = mp.log(2) + g(Delta-1) + g(J+h) - 2*g(Dphi) - 2*g(Dphi-h+1) - g(Delta-h) - g(J+1)
    r += 2*g((J+Delta)/2) + 2*g((Delta-J-d+2)/2) - g(J+Delta-1) - g(Delta-J-d+1)
    r += g((Delta-J-2*d+2*Dphi+2)/2) + g((Delta+J-d+2*Dphi)/2) - g((Delta-J-2*Dphi+2)/2) - g((Delta+J+d-2*Dphi)/2)
    return r

def log_q(Delta, J, Dphi, d=3):
    m2, _ = m2_eta(Delta, J, d)
    return (-mp.log(mp.pi) - log_pMFT(Delta, J, Dphi, d) + mp.log(2*n_J(J, d)) - (d-2)*mp.log(m2)
            + mp.loggamma(Delta-1) + mp.log(2*Delta-d) - mp.loggamma(Delta-d+2))

def log_b(Delta, J, d=3):
    g = mp.loggamma
    h = mp.mpf(d)/2
    return ((d-2)*mp.log(mp.pi) + g(Delta+2-d) + g(J+h-1) + 2*g((J+Delta)/2) + 2*g((Delta+2-d-J)/2)
            - mp.log(2) - g(Delta+1-h) - g(J+d-2) - g(J+Delta) - g(Delta+2-d-J))

def gamma_a(a, nu, d=3):
    h = mp.mpf(d)/2
    return mp.gamma((1+a-h-1j*nu)/2)*mp.gamma((1+a-h+1j*nu)/2)

def a_Dphi(Dphi, d=3):
    h = mp.mpf(d)/2
    return mp.pi**(mp.mpf(d-3)/2)*mp.gamma(Dphi)**2*mp.gamma(Dphi-h+1)**2/(2**(d-5)*mp.gamma(mp.mpf(d-1)/2))

def P(J, eta, d=3):
    """Gegenbauer function normalised to 1 at eta = 1 (CHMRS (2.13)/(B.1)).  J may be complex."""
    return mp.hyp2f1(-J, J+d-2, mp.mpf(d-1)/2, (1-eta)/2)

def dP(J, eta, d=3):
    return mp.diff(lambda x: P(J, x, d), eta)

def Pnu(nu, eta, d=3):
    return P(mp.mpf(2-d)/2 + 1j*nu, eta, d)

def C_series_first(k, nu, Delta, J, Dphi, d=3, printed=False):
    """CHMRS (3.27): (m^{2k}/2) C_{k,nu}[Delta,J] through order 1/m^2.  Returns (leading, leading+first).
    printed=True uses the formula exactly as printed, with (Dphi-d-1)_3; the default uses (2 Dphi-d-1)_3, which is
    what the paper's ancillary file contains (chmrs_series.py) and what the block computation confirms
    (s0_regge_moment.py): the printed Pochhammer symbol is a misprint."""
    m2, eta = m2_eta(Delta, J, d)
    a = 2*Dphi + k - mp.mpf(d)/2
    j = mp.mpf(2-d)/2 + 1j*nu
    P0 = P(j, eta, d); P1 = dP(j, eta, d)
    t1 = P1*((eta**2-1)*(nu**2 - 3*a**2 + 6*a*(2*Dphi-d) - 2) - (d-4)*(d-2))
    t2 = eta*P0*((nu**2 - a*(a+d-2))*(mp.mpf(7*d)/2 + 3*a - 1 - 6*Dphi) - 2*mp.rf((Dphi if printed else 2*Dphi)-d-1, 3)
                 - 2*(a + mp.mpf(d-2)/2)*(1 - a**2 - a*mp.mpf(d-2)/2))
    return P0, P0 + (t1 + t2)/(6*m2)

if __name__ == "__main__":
    mp.mp.dps = 30
    print("# check A: (3.19) -> (3.20), i.e.  4^{2Dphi-1} pi^{d-2} / (b q) * (|p||pbar|)^{-(k+2Dphi-1)}  vs  (2/a) m^{-2k}")
    print("#   with |p||pbar| -> m^2.  ratio -> 1 as m -> infinity at fixed J/Delta")
    for d in (3, 4):
        for Dphi in (mp.mpf('1.2'), mp.mpf('1.7'), mp.mpf(3)):
            row = []
            for Delta, J in ((40, 16), (400, 160), (4000, 1600), (40000, 16000)):
                Delta = mp.mpf(Delta) + mp.mpf('0.3')
                k = 2
                m2, eta = m2_eta(Delta, J, d)
                lhs = mp.exp((2*Dphi-1)*mp.log(4) + (d-2)*mp.log(mp.pi) - log_b(Delta, J, d) - log_q(Delta, J, Dphi, d)
                             - (k+2*Dphi-1)*mp.log(m2))
                rhs = 2/a_Dphi(Dphi, d)/m2**k
                row.append(mp.nstr(lhs/rhs, 10))
            print(f"d={d} Dphi={mp.nstr(Dphi,3)}:", row)
    print("# check B: (3.45) in d = 3:  int_1^inf [d eta] P_j(eta) eta^{-X} = 2^{X+d-4} Gamma((d-1)/2)/(sqrt(pi) Gamma(X)) Gamma((X+j)/2) Gamma((X+2-d-j)/2)")
    d = 3
    for X, nu in ((mp.mpf('3.5'), mp.mpf('0.7')), (mp.mpf(5), mp.mpf('2.3')), (mp.mpf('2.2'), mp.mpf('4.1'))):
        j = mp.mpf(2-d)/2 + 1j*nu
        lhs = mp.quad(lambda e: 2**(d-2)*(e*e-1)**(mp.mpf(d-3)/2)*P(j, e, d)*e**(-X), [1, 2, 10, 100, mp.inf])
        rhs = 2**(X+d-4)*mp.gamma(mp.mpf(d-1)/2)/(mp.sqrt(mp.pi)*mp.gamma(X))*mp.gamma((X+j)/2)*mp.gamma((X+2-d-j)/2)
        print(f"X={mp.nstr(X,3)} nu={mp.nstr(nu,3)}: lhs={mp.nstr(lhs,15)} rhs={mp.nstr(rhs,15)}")
