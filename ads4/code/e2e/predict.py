"""
Predicted leading heavy-block value of the tensor Regge moments (derivation: SPINNING_E2E.md sect. 4).

   (-1)^J Pihat^{(q1,q3)}_{k,nu}[ stripped W_diff ]  ->
        (pi/4) 4^{l+1+Dphi} / b^{(1-Dphi)}_{Delta,J} * gamma_{a_x}(nu) gamma_{a_y}(nu) m^{-(a_x+a_y)} * T_l * sum_{h1,h3} M^{q1q3}_{h1h3}(nu; a_y) F^nu_{mu}(beta),
   a_x = 2 Dphi + k - 1,  a_y = 2(l+1) + k - 1,
   T_1 = -(1/2) (m^3 sinh(beta)/2)^2,     T_2 = (m^6 sinh^2(beta)/8)^2,
   M = R + g1 S (laplace.py),  g1 = gamma_{a_y - 1}(nu)/gamma_{a_y}(nu),  mu = harmonic returned by laplace.M_coeff
   (0 for helicity preserved, +-2l for helicity flip).
W_diff is the partial wave built with the conservation coefficients c_w of structures.py; in terms of the standard conserved
structure (S_0 coefficient 1) the left/right structures are 2^{-J} L_0 times the standard one.
"""
import mpmath as mp, sympy as sp
from formulas import gamma_a, log_b, Fnu_mp, m2_eta
from laplace import M_table, nu as symnu, a as syma, g1 as symg1

_Mnum = {}
def Mfun(l):
    if l not in _Mnum:
        T = M_table(l)
        _Mnum[l] = {k: (sp.lambdify((symnu, syma), v[0], 'mpmath'), sp.lambdify((symnu, syma), v[1], 'mpmath'), v[2]) for k, v in T.items()}
    return _Mnum[l]

def components(l, k, nu, q1, q3):
    """(M_pres, M_flip): coefficients of F_0 = P and of F_{2l} in sum_{h1,h3} M F  (a = a_y = 2l + k + 1)."""
    ay = 2*(l+1) + k - 1
    g1v = gamma_a(ay - 1, nu)/gamma_a(ay, nu)
    T = Mfun(l)
    pres = flip = mp.mpc(0)
    for h1 in (1, -1):
        for h3 in (1, -1):
            R, S, mu = T[(q1, q3, h1, h3)]
            val = R(nu, ay) + g1v*S(nu, ay)
            if mu == 0: pres += val
            else:
                assert abs(mu) == 2*l, mu
                flip += val
    return pres, flip

def scalar_prefactor(l, Delta, J, Dphi, k, nu):
    Delta = mp.mpf(Delta)
    m2, eta = m2_eta(Delta, J)
    beta = mp.acosh(eta)
    ax, ay = 2*Dphi + k - 1, 2*(l+1) + k - 1
    lg = mp.log(mp.pi/4) + (l + 1 + Dphi)*mp.log(4) - log_b(Delta, J, 1 - Dphi) - (ax + ay)/2*mp.log(m2)
    m = mp.sqrt(m2)
    if l == 1:
        T = -(m**3*mp.sinh(beta)/2)**2/2
    else:
        T = (m**6*mp.sinh(beta)**2/8)**2
    return mp.exp(lg)*gamma_a(ax, nu)*gamma_a(ay, nu)*T, beta

def prediction(l, Delta, J, Dphi, k, nu, q1, q3, parts=False):
    pref, beta = scalar_prefactor(l, Delta, J, Dphi, k, nu)
    pres, flip = components(l, k, nu, q1, q3)
    P0 = Fnu_mp(nu, 0, beta); Pf = Fnu_mp(nu, 2*l, beta)
    if parts:
        return pref*pres*P0, pref*flip*Pf
    return pref*(pres*P0 + flip*Pf)
