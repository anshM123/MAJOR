"""
Scalar-skeleton AdS remainder relative to the hybrid (flat + leading-Regge curvature) description, per heavy state (Delta, J)
and per transverse momentum nu = Lam p (Df = 3: the <TTTT> skeleton, a = k + 5).

Definitions (m^{2k}-normalised heavy densities, x = nu/m, eta = eta_AdS, JJ = J(J+1)):
  E_k(x)    = (2 - x^2)/(1 - x^2)^{1+k/2}                     scalar bulk-point envelope, BP_k = E_k(x) P_J(1 - 2x^2)  [CHMRS (3.40)]
  hyb_k     = BP_k + 2[P_{-1/2+i nu}(eta) - J_0((2J+1) x)]       (flat + [Regge - its flat limit at fixed bhat])
  exact_k   = BP_k + r1(x,J)/m^2 + r2(x,J)/m^4 + ...             fixed-J series of CHMRS (ancillary, order 8)
  residual  rho_k = (exact_k - hyb_k)/E_k(x)  =  [r1res/m^2 + r2res/m^4]/E_k,   r_ires = r_i - (Regge-curvature part, exact)
 r1res, r2res are truncated power series in x^2 with coefficients polynomial in JJ (cache/residual_coeffs_Df3.pkl, t1).
 'controlled' at a node if the truncated series is decreasing (|t_last| <= 0.8 |t_prev|); error = geometric tail.
 Fallback for nodes not controlled (Regge regime, x small): the fixed-eta series seriesC (order 8, validated in v1):
  rho_k = (S_k - hyb_k)/E_k, controlled if its last term is decreasing and < 2 % of E_k.
Forward terms at nu = 0 (fixed-eta series, uniform in beta):
  R4 = -m^2 d/dnu^2 [m^8 C_4]|_0 - (2 m^2 W(beta) - 11/2),   R5 = m^10 C_{5,0} - 2 P_{-1/2}(cosh beta)
 (the hybrid values 2m^2W - 11/2 and 2P0 are flat (2JJ - 5, 2) + leading-Regge curvature).
"""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, mpmath as mp
from scipy.special import j0, eval_legendre
import skelfast

DF = float(os.environ.get('UNIFORM_DF', '3'))       # skeleton external dimension (3: leading-order <TTTT> choice; 2: SPINNING_MFT choice)
TAG = '' if DF == 3 else '_Df%g' % DF
RC = pickle.load(open(os.path.join(HERE, 'cache', 'residual_coeffs_Df%g.pkl' % DF), 'rb'))

def Ek(k, x):
    return (2 - x*x)/(1 - x*x)**(1 + k/2)

def Ek_tail(k, x, jmax=4):
    """E_k(x) minus its Taylor polynomial in y = x^2 through y^jmax (>= 0)."""
    from math import factorial
    y = x*x
    C = lambda j: 0.0 if j < 0 else np.prod([1 + k/2 + i for i in range(j)])/factorial(j)
    poly = sum((2*C(j) - C(j - 1))*y**j for j in range(jmax + 1))
    return np.abs(Ek(k, x) - poly)

def _series_terms(coefs, x2, JJ, offset):
    """terms t_j = x2^{j-offset} * poly_j(JJ)"""
    out = []
    for i, c in enumerate(coefs):
        out.append(x2**i*np.polynomial.polynomial.polyval(JJ, c))
    return np.array(out)

def _tail(t):
    """sum, error, controlled flag for a truncated series with terms t (rows) at each node."""
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
    t1 = _series_terms(RC[k]['r1'], x2, JJ, 0)
    t2 = _series_terms(RC[k]['r2'], x2, JJ, 0)
    s1, e1, c1, l1 = _tail(t1)
    s2, e2, c2, l2 = _tail(t2)
    E = Ek(k, x)
    # a last term that is negligible in absolute terms also counts as controlled (accidental near-zeros of t_{last-1})
    neg1 = l1/m2 <= 1e-4*E; neg2 = l2/m2**2 <= 1e-4*E
    e1 = np.where(c1, e1, 2*l1); e2 = np.where(c2, e2, 2*l2)
    c1 = c1 | neg1; c2 = c2 | neg2
    rho = (s1/m2 + s2/m2**2)/E
    # 1/m^6 budget (estimate, after the referee): next order ~ (s2/m^4) * (|s2|/max(|s1|,1))/m^2
    e6 = np.abs(s2)/m2**2*(np.abs(s2)/np.maximum(np.abs(s1), 1.0))/m2
    err = (e1/m2 + e2/m2**2 + e6)/E
    # require the 1/m^4 piece not to dominate the 1/m^2 piece by more than the expected 1/m^2 suppression
    ctrl = c1 & c2 & (np.abs(s2/m2**2) <= np.maximum(np.abs(s1/m2), 1e-3*E))
    return rho, err, ctrl

def hyb_regge_part(nus, eta, J, m):
    """2[P_{-1/2+i nu}(eta) - J0((2J+1) nu/m)] (exact, mpmath P)."""
    P = skelfast.Pders_vec(nus, eta, nmax=0)[0]
    return 2*(P - j0((2*J + 1)*nus/m))

def residual_seriesC(k, nus, Delta, J, m2, eta, Pd):
    T = skelfast.series_terms_vec(k, DF, nus, eta, Pd=Pd)
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
    # error: tail of the 1/m^2 series + truncation of the bulk-point x-series beyond x^8 (bounded by the x^10 term of E_k)
    # truncation of the bulk-point x-series beyond x^8 inside the order-8 fixed-eta series: twice the exact tail of E_k
    err = (last*q/(1 - np.minimum(q, 0.8)) + 2*Ek_tail(k, x))/E
    return rho, err, ctrl

def forward_remainders(Delta, J, m2, eta, W, P0, h=0.05):
    """R4, R5 and error flags from the fixed-eta series at nu = 0."""
    nus = np.array([0.0, h, np.sqrt(2.0)*h])      # nu^2 = 0, h^2, 2h^2 (equally spaced in nu^2)
    Pd = skelfast.Pders_vec(nus, eta)
    T4 = skelfast.series_terms_vec(4, DF, nus, eta, Pd=Pd)
    T5 = skelfast.series_terms_vec(5, DF, nus[:1], eta, Pd=Pd[:, :1])
    tm4 = np.array([T4[j]/m2**j for j in range(5)])           # (5, 3)
    d4 = (-3*tm4[:, 0] + 4*tm4[:, 1] - tm4[:, 2])/(2*h*h)       # d/dnu^2 per order
    terms4 = -m2*d4
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
    return R4, e4, c4, R5, e5, c5, terms4, tm5
