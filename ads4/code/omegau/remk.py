"""Scalar-skeleton remainder rho_k (k = 2..7) relative to the per-k hybrid (bulk point + leading-Regge curvature), per heavy state and node.

Adapted copy of uniform/remainder.py + uniform/skelfast.py (same definitions, same heuristic control cut [M]), extended to k = 2..7,
with caches in omegau/cache/.  Skeleton external dimension from env OMEGAU_DF (3: leading-order <TTTT> choice, 2: MFT choice, 2.5).
  rho_k = (exact_k - hyb_k)/E_k,  hyb_k = E_k(x) P_J(1-2x^2) + 2[P_{-1/2+i nu}(eta) - J_0((2J+1)x)],  E_k = (2-x^2)/(1-x^2)^{1+k/2}
  exact_k: (1) fixed-J series BP + r1/m^2 + r2/m^4 (curvature-subtracted coefficients of u2), controlled if each truncated x^2-series has
           last/previous ratio <= 0.8 and the 1/m^4 piece does not dominate; error = geometric tail + 1/m^6 budget (estimate);
           (2) else, if x <= 0.6: fixed-eta series (CHMRS seriesC, order 1/m^8), controlled if decreasing and last term < 2% E_k;
               error = tail + 2 x (exact tail of E_k beyond x^8);   (3) else uncontrolled.
Forward (nu = 0): R4 = -m^2 d_{nu^2}(m^8 C_4)|_0 - (2 m^2 W - 11/2),  R5 = m^10 C_{5,0} - 2 P_{-1/2}(cosh beta).
"""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
import numpy as np, mpmath as mp
from scipy.special import j0, eval_legendre

DF = float(os.environ.get('OMEGAU_DF', '3'))
TAG = 'Df%g' % DF
RC = pickle.load(open(os.path.join(HERE, 'cache', 'rcoef_%s.pkl' % TAG), 'rb'))
KS = (2, 3, 4, 5, 6, 7)

def Ek(k, x):
    return (2 - x*x)/(1 - x*x)**(1 + k/2)

def Ek_tail(k, x, jmax=4):
    from math import factorial
    y = x*x
    C = lambda j: 0.0 if j < 0 else np.prod([1 + k/2 + i for i in range(j)])/factorial(j)
    poly = sum((2*C(j) - C(j - 1))*y**j for j in range(jmax + 1))
    return np.abs(Ek(k, x) - poly)

# ---------------------------------------------------------------- fixed-eta series (copy of skelfast with own cache)
def _build(k_, Df_):
    import sympy as sp
    from chmrs_series import load, d, k, Df, al, nu, eta, P     # reads uniform/cache/chmrs_series.pkl (read-only)
    s = load()
    Df_s = sp.nsimplify(Df_)
    arrs = []
    for c in s['seriesC']:
        e = sp.expand(c.subs({d: 3, k: k_, Df: Df_s, al: 2*Df_s + k_ - 1}))
        per_n = []
        for n in range(5):
            cn = sp.expand(e.coeff(P[n]))
            if cn == 0:
                per_n.append(None); continue
            nu2 = sp.Symbol('nu2')
            cn = sp.together(sp.expand(cn.subs(nu, sp.sqrt(nu2))))
            num, den = sp.fraction(cn)
            mats = []
            for q in (sp.expand(num), sp.expand(den)):
                pol = sp.Poly(q, nu2, eta)
                A = np.zeros((pol.degree(nu2) + 1, pol.degree(eta) + 1))
                for (i, j), v in pol.terms():
                    A[i, j] = float(v)
                mats.append(A)
            per_n.append(tuple(mats))
        arrs.append(per_n)
    return arrs

_C = {}
def coeffs(k, Df=DF):
    key = (k, float(Df))
    if key not in _C:
        pk = os.path.join(HERE, 'cache', 'sf_%d_%g.pkl' % key)
        if os.path.exists(pk):
            _C[key] = pickle.load(open(pk, 'rb'))
        else:
            _C[key] = _build(k, Df)
            pickle.dump(_C[key], open(pk, 'wb'))
    return _C[key]

def Pders_vec(nus, eta, nmax=4):
    out = np.empty((nmax + 1, len(nus)))
    z = (1 - mp.mpf(eta))/2
    for i, nu in enumerate(nus):
        a, b = mp.mpf(1)/2 + 1j*mp.mpf(nu), mp.mpf(1)/2 - 1j*mp.mpf(nu)
        for n in range(nmax + 1):
            pre = (-mp.mpf(1)/2)**n*mp.rf(a, n)*mp.rf(b, n)/mp.factorial(n)
            out[n, i] = float((pre*mp.hyp2f1(a + n, b + n, 1 + n, z)).real)
    return out

def series_terms_vec(k, nus, eta, Pd):
    nus = np.asarray(nus, float)
    nu2 = nus**2
    T = []
    for per_n in coeffs(k):
        t = np.zeros_like(nus)
        for n, AB in enumerate(per_n):
            if AB is None:
                continue
            A, B = AB
            ev = lambda M: np.polynomial.polynomial.polyval(nu2, np.polynomial.polynomial.polyval(eta, M.T))
            t += Pd[n]*ev(A)/ev(B)
        T.append(t)
    return np.array(T)

# ---------------------------------------------------------------- remainders
def _series_terms(coefs, x2, JJ):
    return np.array([x2**i*np.polynomial.polynomial.polyval(JJ, c) for i, c in enumerate(coefs)])

def _tail(t):
    s = t.sum(axis=0)
    last, prev = np.abs(t[-1]), np.abs(t[-2])
    with np.errstate(divide='ignore', invalid='ignore'):
        q = np.where(prev > 0, last/prev, np.where(last > 0, np.inf, 0.0))
    ctrl = q <= 0.8
    err = np.where(ctrl, last*q/(1 - np.minimum(q, 0.8) + 1e-300), np.inf)
    return s, err, ctrl, last

def residual_fixedJ(k, x, J, m2):
    JJ = J*(J + 1.0)
    x2 = x*x
    s1, e1, c1, l1 = _tail(_series_terms(RC[k]['r1'], x2, JJ))
    s2, e2, c2, l2 = _tail(_series_terms(RC[k]['r2'], x2, JJ))
    E = Ek(k, x)
    neg1 = l1/m2 <= 1e-4*E; neg2 = l2/m2**2 <= 1e-4*E
    e1 = np.where(c1, e1, 2*l1); e2 = np.where(c2, e2, 2*l2)
    c1 = c1 | neg1; c2 = c2 | neg2
    rho = (s1/m2 + s2/m2**2)/E
    e6 = np.abs(s2)/m2**2*(np.abs(s2)/np.maximum(np.abs(s1), 1.0))/m2
    err = (e1/m2 + e2/m2**2 + e6)/E
    ctrl = c1 & c2 & (np.abs(s2/m2**2) <= np.maximum(np.abs(s1/m2), 1e-3*E))
    return rho, err, ctrl

def residual_seriesC_eta(k, nus, J, m2, eta, Pd):
    T = series_terms_vec(k, nus, eta, Pd)
    terms = np.array([T[j]/m2**j for j in range(5)])
    S = terms.sum(axis=0)
    m = np.sqrt(m2); x = nus/m
    BP = Ek(k, x)*eval_legendre(J, 1 - 2*x*x)
    hyb = BP + 2*(Pd[0] - j0((2*J + 1)*x))
    E = Ek(k, x)
    last, prev = np.abs(terms[-1]), np.abs(terms[-2])
    q = np.where(prev > 0, last/np.maximum(prev, 1e-300), 1.0)
    ctrl = ((q <= 0.8) & (last <= 0.02*E)) | (last <= 1e-5*E)
    q = np.where(q <= 0.8, q, 0.5)
    rho = (S - hyb)/E
    err = (last*q/(1 - np.minimum(q, 0.8)) + 2*Ek_tail(k, x))/E
    return rho, err, ctrl

def forward_remainders(J, m2, eta, W, P0, h=None):
    if h is None:      # nu-step on the scale 1/beta of the Regge kernel (avoids round-off of m^2 d/dnu^2 at large m); 0.05 for beta >= 2
        beta = float(np.arccosh(eta)) if eta > 1 else 1e-12
        h = min(max(0.05, 0.1/beta), 0.02*np.sqrt(m2))
    nus = np.array([0.0, h, np.sqrt(2.0)*h])
    Pd = Pders_vec(nus, eta)
    T4 = series_terms_vec(4, nus, eta, Pd)
    T5 = series_terms_vec(5, nus[:1], eta, Pd[:, :1])
    tm4 = np.array([T4[j]/m2**j for j in range(5)])
    terms4 = -m2*(-3*tm4[:, 0] + 4*tm4[:, 1] - tm4[:, 2])/(2*h*h)
    S4 = terms4.sum()
    tm5 = np.array([T5[j][0]/m2**j for j in range(5)])
    S5 = tm5.sum()
    R4 = S4 - (2*m2*W - 5.5)
    R5 = S5 - 2*P0
    def ctl(t, scale):
        q = abs(t[-1])/max(abs(t[-2]), 1e-300)
        return (q <= 0.8) and abs(t[-1]) <= 0.05*scale, abs(t[-1])*q/(1 - min(q, 0.8))
    c4, e4 = ctl(terms4, max(abs(S4), 1.0))
    c5, e5 = ctl(tm5, 2.0)
    return R4, e4, c4, R5, e5, c5

def skeleton_state(st, J, x, W, P0):
    """per-k (rho, err, ctrl) at nodes x = nu/m, and forward (R4, e4, R5, e5, ctrl) for the heavy state st (dict from common.state)."""
    m2, eta, m = st['m2'], st['eta'], st['m']
    nus = x*m
    out = {}
    Pd = None
    for k in KS:
        rho, err, ctrl = residual_fixedJ(k, x, J, m2)
        need = ~ctrl & (x <= 0.6)
        if need.any():
            if Pd is None:
                Pd = np.zeros((5, len(x)))
                idx = np.where(x <= 0.6)[0]
                Pd[:, idx] = Pders_vec(nus[idx], eta)
            r2, e2, c2 = residual_seriesC_eta(k, nus, J, m2, eta, Pd)
            use = need & c2
            rho = np.where(use, r2, rho); err = np.where(use, e2, err); ctrl = ctrl | use
        out[k] = (rho, np.where(ctrl, err, 0.0), ctrl)
    R4, e4, c4, R5, e5, c5 = forward_remainders(J, m2, eta, W, P0)
    out['fw'] = (R4, e4, R5, e5, bool(c4 and c5))
    return out
