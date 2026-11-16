"""
Predicted leading heavy value of the doubly-spinning tensor Regge moments (factorised formula; derivation SPINNING_4PT.md sect. 3):

   (-1)^J Pihat^q_{k,nu}[ s13^{l+1} s24^{l+1} W[v_a, v_b] ]  ->
        (pi/4) 4^{2(l+1)} / b^{(0)}_{Delta,J} * gamma_a(nu)^2 m^{-2a} * sum_h T^{ab}_h M^{q1q3}_{h1h3}(nu; a) M^{q2q4}_{h2h4}(-nu; a) K^nu_{muL(h), muR(h)}(beta),
   a = 2l + k + 1,  T^{ab}_h = V_a(h1, h2) V_b(h3, h4)  (heavy vertex tensors, vertex2 / t3_vertex_l1),
   M = R + g1 S from e2e/laplace.py (one Regge pair each; the x-pair with nu -> -nu),  K = moments4.kernel (two-pair harmonic function).
Heavy vertex tensors (l = 1; A_+ = Delta - 3/2, A_- = J + 1/2, m^2 = A_+^2 - A_-^2), in the frames of moments4 transported to
(pbarhat, phat):   V_1(h1,h2) = (A_+ A_-)^2 m^2/2 * [1 (h2 = -h1: 'same' helicity) | 1 (h2 = h1: 'opposite')],
                   V_2(h1,h2) = -(1/2) * [A_+^2 (same) | -A_-^2 (opposite)].
"""
import os, sys, functools
import mpmath as mp, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
from formulas import gamma_a, log_b                      # noqa: E402
from laplace import M_table, nu as symnu, a as syma      # noqa: E402
from moments4 import kernel_mp

@functools.lru_cache(None)
def Mfun(l):
    T = M_table(l)
    return {k: (sp.lambdify((symnu, syma), v[0], 'mpmath'), sp.lambdify((symnu, syma), v[1], 'mpmath'), v[2]) for k, v in T.items()}

def Mval(l, q1, q3, h1, h3, nu, a):
    R, S, mu = Mfun(l)[(q1, q3, h1, h3)]
    g1 = gamma_a(a - 1, nu)/gamma_a(a, nu)
    return R(nu, a) + g1*S(nu, a), mu

def weyl(Delta, J):
    D = mp.mpf(Delta)
    Ap, Am = D - mp.mpf(3)/2, J + mp.mpf(1)/2
    m = mp.sqrt(Ap*Ap - Am*Am)
    beta = 2*mp.atanh(Am/Ap)
    return Ap, Am, m, beta

def V_l1(a, Ap, Am, m, h1, h2):
    same = (h2 == -h1)
    if a == 1:
        return (Ap*Am)**2*m**2/2
    if a == 2:
        return -(Ap**2 if same else -Am**2)/2
    raise ValueError

def T_l1(cL, cR, Ap, Am, m, h):
    """T_h for coefficient vectors cL, cR on (v_1, v_2) (numbers)."""
    h1, h2, h3, h4 = h
    VL = sum(cL[i]*V_l1(i+1, Ap, Am, m, h1, h2) for i in range(2))
    VR = sum(cR[i]*V_l1(i+1, Ap, Am, m, h3, h4) for i in range(2))
    return VL*VR

_Kc = {}
def Kcached(nu, muL, muR, beta):
    key = (float(nu), muL, muR, mp.nstr(beta, 25))
    if key not in _Kc:
        _Kc[key] = kernel_mp(nu, muL, muR, beta)
    return _Kc[key]

def prefactor(l, Delta, J, k, nu):
    Ap, Am, m, beta = weyl(Delta, J)
    a = 2*l + k + 1
    lg = mp.log(mp.pi/4) + 2*(l+1)*mp.log(4) - log_b(mp.mpf(Delta), J, 0) - 2*a*mp.log(m)
    return mp.exp(lg)*gamma_a(a, nu)**2

def prediction(l, Delta, J, k, nu, q, cL, cR, parts=False, Tfun=None):
    """(-1)^J Pihat^q predicted.  q = (q1, q3, q2, q4).  cL, cR: coefficients on (v_1, v_2).  parts: dict (muL, muR) -> contribution."""
    assert l == 1 or Tfun is not None
    Ap, Am, m, beta = weyl(Delta, J)
    a = 2*l + k + 1
    q1, q3, q2, q4 = q
    pre = prefactor(l, Delta, J, k, nu)
    acc = {}
    for h1 in (1, -1):
        for h3 in (1, -1):
            ML, muL = Mval(l, q1, q3, h1, h3, nu, a)
            if ML == 0: continue
            for h2 in (1, -1):
                for h4 in (1, -1):
                    MR, muR = Mval(l, q2, q4, h2, h4, -nu, a)
                    if MR == 0: continue
                    T = (Tfun or T_l1)(cL, cR, Ap, Am, m, (h1, h2, h3, h4))
                    acc[(muL, muR)] = acc.get((muL, muR), 0) + T*ML*MR
    out = {key: pre*v*Kcached(nu, key[0], key[1], beta) for key, v in acc.items()}
    if parts:
        return out
    return sum(out.values())
