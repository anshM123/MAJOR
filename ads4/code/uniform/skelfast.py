"""
Vectorised evaluation of the CHMRS order-8 series m^{2k} C_{k,nu}[Delta,J] (seriesC) for d = 3 and fixed (k, Df):
each order-j coefficient is  sum_n P^{(n)}(eta) * poly_{j,n}(nu^2, eta)  with exact rational polynomial coefficients
(precomputed with sympy from the ancillary file; cached in cache/skelfast_<k>_<Df>.pkl).
"""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, mpmath as mp

def _build(k_, Df_):
    import sympy as sp
    from chmrs_series import load, d, k, Df, al, nu, eta, P
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
def coeffs(k, Df):
    key = (k, float(Df))
    if key not in _C:
        pk = os.path.join(HERE, 'cache', 'skelfast_%d_%g.pkl' % key)
        if os.path.exists(pk):
            _C[key] = pickle.load(open(pk, 'rb'))
        else:
            _C[key] = _build(k, Df)
            pickle.dump(_C[key], open(pk, 'wb'))
    return _C[key]

def Pders_vec(nus, eta, nmax=4):
    """array (nmax+1, len(nus)) of eta-derivatives of P_{-1/2+i nu}(eta)."""
    out = np.empty((nmax + 1, len(nus)))
    z = (1 - mp.mpf(eta))/2
    for i, nu in enumerate(nus):
        a, b = mp.mpf(1)/2 + 1j*mp.mpf(nu), mp.mpf(1)/2 - 1j*mp.mpf(nu)
        for n in range(nmax + 1):
            pre = (-mp.mpf(1)/2)**n*mp.rf(a, n)*mp.rf(b, n)/mp.factorial(n)
            out[n, i] = float((pre*mp.hyp2f1(a + n, b + n, 1 + n, z)).real)
    return out

def series_terms_vec(k, Df, nus, eta, Pd=None):
    """terms t_j (j = 0..4) of m^{2k} C_{k,nu} = sum_j t_j / m^{2j} at eta_AdS = eta; arrays over nus."""
    nus = np.asarray(nus, float)
    if Pd is None:
        Pd = Pders_vec(nus, eta)
    nu2 = nus**2
    T = []
    for per_n in coeffs(k, Df):
        t = np.zeros_like(nus)
        for n, AB in enumerate(per_n):
            if AB is None:
                continue
            A, B = AB      # numerator / denominator: polys in nu2 (rows) and eta (cols)
            ev = lambda M: np.polynomial.polynomial.polyval(nu2, np.polynomial.polynomial.polyval(eta, M.T))
            t += Pd[n]*ev(A)/ev(B)
        T.append(t)
    return np.array(T)
