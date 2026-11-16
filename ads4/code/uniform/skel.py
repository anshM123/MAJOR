"""
Scalar-skeleton heavy moments, d = 3 (CHMRS = arXiv:2106.10274):

 * exact asymptotic (Regge) series of CHMRS's ancillary file to order 1/m^8:  m^{2k} C_{k,nu}[Delta,J]  (seriesC) and
   a_{Df} m^{2k}/gamma^2 Pihat_{k,nu}[Delta,J] (seriesPi), as numeric functions;  al = 2 Df + k - 1 (identified in v1);
 * the fixed-J (bulk-point) series seriesCfixedJ;
 * the leading bulk-point formula  C^BP_{k,nu} = (2m^2-nu^2)/(m^2-nu^2) P_J(1-2nu^2/m^2)/(m^2(m^2-nu^2))^{k/2}  [CHMRS (3.40)];
 * the large-(Delta,J,nu) saddle formula  S_nu  [CHMRS app. E, (E.2)-(E.3)];
 * normalisation factors q_{Delta,J}, a_{Df}, gamma_a (CHMRS (2.43)-(2.44), (2.55), (2.57)).
"""
import os, sys, pickle, functools
import numpy as np, mpmath as mp, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# ---------------------------------------------------------------- kinematics
def m2_eta(Delta, J):
    a, ab = Delta - J - 2, Delta + J - 1
    return a*ab, (a*a + ab*ab)/(2*a*ab)

def Delta_of(m2, J):
    """Delta with given m^2 = (Delta-J-2)(Delta+J-1) and J."""
    # (D - J - 2)(D + J - 1) = m2 ;  D^2 - 3D + (2 - J - J^2 + J... ) : solve
    A, B = -J - 2, J - 1        # (D + A)(D + B) = m2
    s = (A + B)
    return (-s + np.sqrt(s*s - 4*(A*B - m2)))/2

# ---------------------------------------------------------------- harmonic function and derivatives
def Pder(nu, eta, n):
    """n-th eta-derivative of P_{-1/2+i nu}(eta) (real)."""
    nu = mp.mpf(nu); eta = mp.mpf(eta)
    a, b = mp.mpf(1)/2 + 1j*nu, mp.mpf(1)/2 - 1j*nu
    pre = (-mp.mpf(1)/2)**n*mp.rf(a, n)*mp.rf(b, n)/mp.factorial(n)
    return (pre*mp.hyp2f1(a + n, b + n, 1 + n, (1 - eta)/2)).real

def Pders(nu, eta, nmax=4):
    return [float(Pder(nu, eta, n)) for n in range(nmax + 1)]

# ---------------------------------------------------------------- normalisations (d = 3)
def gamma_a(a, nu):
    return abs(mp.gamma((a - mp.mpf(1)/2 + 1j*nu)/2))**2

def a_Df(Df):
    Df = mp.mpf(Df)
    return mp.gamma(Df)**2*mp.gamma(Df - mp.mpf(1)/2)**2/(mp.mpf(2)**(-2)*mp.gamma(1))

def log_pMFT(Delta, J, Df):
    g = mp.loggamma; D = mp.mpf(Delta); Df = mp.mpf(Df); d = 3
    return (mp.log(2) + g(D - 1) + g(J + mp.mpf(d)/2) - 2*g(Df) - 2*g(Df - mp.mpf(d - 2)/2) - g(D - mp.mpf(d)/2) - g(J + 1)
            + 2*g((J + D)/2) + 2*g((D - J - d + 2)/2) - g(J + D - 1) - g(D - J - d + 1)
            + g((D - J - 2*d + 2*Df + 2)/2) + g((D + J - d + 2*Df)/2) - g((D - J - 2*Df + 2)/2) - g((D + J + d - 2*Df)/2))

def log_q(Delta, J, Df):
    """q_{Delta,J} of CHMRS (2.43), d = 3, D = 4: n_J = 16 pi (2J+1)."""
    D = mp.mpf(Delta); d = 3
    m2 = (D - J - d + 1)*(D + J - 1)
    nJ = 16*mp.pi*(2*J + 1)
    return (-mp.log(mp.pi) - log_pMFT(Delta, J, Df) + mp.log(2*nJ) - (d - 2)*mp.log(m2)
            + mp.loggamma(D - 1) + mp.log(2*D - d) - mp.loggamma(D - d + 2))

def log_b(Delta, J):
    g = mp.loggamma; D = mp.mpf(Delta)
    return (mp.log(mp.pi) + g(D - 1) + g(J + mp.mpf(1)/2) + 2*g((J + D)/2) + 2*g((D - 1 - J)/2)
            - mp.log(2) - g(D - mp.mpf(1)/2) - g(J + 1) - g(J + D) - g(D - 1 - J))

# ---------------------------------------------------------------- CHMRS series, numeric
_SER = None
def _series():
    global _SER
    if _SER is None:
        pk = os.path.join(HERE, 'cache', 'series_num.pkl')
        if os.path.exists(pk):
            _SER = pickle.load(open(pk, 'rb'))
        else:
            from chmrs_series import load, d, k, Df, al, nu, eta, P, J2
            s = load()
            out = {}
            for name in ('seriesPi', 'seriesC', 'seriesCfixedJ'):
                cs = []
                for c in s[name]:
                    c3 = sp.expand(c.subs(d, 3))
                    cs.append(sp.srepr(c3))
                out[name] = cs
            _SER = out
            pickle.dump(out, open(pk, 'wb'))
    return _SER

@functools.lru_cache(None)
def _lamb(name, j):
    from chmrs_series import k, Df, al, nu, eta, P, J2
    expr = sp.sympify(eval(_series()[name][j], sp.__dict__.copy()))  # srepr round trip
    args = [k, Df, al, nu, eta, J2] + list(P[:5])
    return sp.lambdify(args, expr, 'math')

def series_terms(name, k, Df, nu, eta, J=0, Pd=None):
    """list of the 5 terms (order m2^0 .. m2^-4 coefficients) of seriesPi/seriesC at given (k, Df, nu, eta);
       J only for seriesCfixedJ. al = 2 Df + k - 1."""
    al = 2*Df + k - 1
    if Pd is None and name != 'seriesCfixedJ':
        Pd = Pders(nu, eta)
    if Pd is None:
        Pd = [0.0]*5
    J2 = J*(J + 1)
    return [_lamb(name, j)(k, Df, al, nu, eta, J2, *Pd) for j in range(5)]

def series_value(name, k, Df, nu, Delta, J, order=4, Pd=None):
    m2, eta = m2_eta(Delta, J)
    t = series_terms(name, k, Df, nu, eta, J, Pd)
    return sum(t[j]/m2**j for j in range(order + 1))

# ---------------------------------------------------------------- leading closed forms
def legP(J, x):
    return float(mp.legendre(J, x))

def C_BP(k, nu, m2, J):
    """CHMRS (3.40)/(3.39): m^{2k} C^BP = (2-x^2)/(1-x^2)^{1+k/2} P_J(1-2x^2), x = nu/m (returned without the m^{-2k})."""
    x2 = nu*nu/m2
    return (2 - x2)/(1 - x2)**(1 + k/2)*legP(J, 1 - 2*x2)

def C_regge_lead(nu, eta):
    """m^{2k} C^Regge = 2 P_{-1/2+i nu}(eta_AdS)."""
    return 2*float(Pder(nu, eta, 0))

def S_formula(k, nu, Delta, J, d=3):
    """m^{2k} x CHMRS (E.2): (-1)^J (S_nu + c.c.) with S_nu of (E.3); tau = Delta-J-d+1, taub = Delta+J-1."""
    t, tb, nu = mp.mpf(Delta - J - d + 1), mp.mpf(Delta + J - 1), mp.mpf(nu)
    st, stb = mp.sqrt(t*t - nu*nu), mp.sqrt(tb*tb - nu*nu)
    I = 1j
    pre = 2**(d - 3)*mp.gamma(mp.mpf(d - 1)/2)/mp.sqrt(mp.pi)
    f1 = (I*t*t*tb*tb/(nu*(tb*tb - t*t)))**(mp.mpf(d - 2)/2)
    f2 = (t*tb + st*stb)/((t*tb)**(mp.mpf(k + 1)/2)*(t*t - nu*nu)**(mp.mpf(k + 3)/4)*(tb*tb - nu*nu)**(mp.mpf(k + 3)/4))
    f3 = ((st + I*nu)/t)**t*((stb - I*nu)/tb)**tb*((st + t)/(stb + tb))**(I*nu)
    f4 = ((tb*st + t*stb)/2)**(mp.mpf(d)/2 - 2)
    S = pre*f1*f2*f3/f4
    m2 = t*tb
    return float(((-1)**J*(S + mp.conj(S))).real*m2**k)
