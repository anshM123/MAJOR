"""p4win/p4w_contour.py (P4_WINDOW.md §2.3): constants of the complement (contour) part of Lemma P4W.

For a correction G on [p_w, 1] (p = nu/Lam) of class AG, G = G_+ + G_-, |G_+-(p +- i y)| <= M(Re p) exp(-kappa(Re p) y), the k-th moment
cost 2 int_{p_w}^1 psi_k G dp (k = 2; k = 3 with weight u sigma, no factor 2) obeys, by Cauchy on the rectangles [p_w, 1] x [0, +-T]:
    |2 int psi G| <= 2 * 2 * [ V(p_w) M(p_w) + V(1) M(1) + Top ],
    V(p0) = int_0^T |psi(p0 + i y)| exp(-kappa(p0) y) dy,   Top = exp(-kappa_min T) int_{p_w}^1 |psi(p + i T)| M(p) dp.
(the first factor 2 is the ring weight of the k = 2 action, the second the two waves).  kappa: 'uniform' = Lam beta_*, beta_* = 0.97 beta
(HLJ_ANALYTIC §3: AG to all orders with beta_* >= 1.002 beta); 'local' = Lam beta_*/sqrt(1 - x^2) (the local frequency of the flat
partial-wave phase, as HLJ's E_edge constant) on 0 <= y <= T1F/kappa (T1F = 15 default; env P4W_T1), then the uniform rate (continuous majorant); rectangle height
T = 25/(0.97 Lam beta).  psi_k are evaluated as complex polynomials (exact NF02 coefficients).
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import p4w_common as C

F = C.NF
T1F = float(os.environ.get('P4W_T1', '15'))   # local rate kappa(p) on 0 <= y <= T1F/kappa(p), uniform rate beyond (AG' premise)


def psi_c(k, p):
    p = np.asarray(p, complex)
    co = F.X if k == 2 else F.Y
    return (1 - p)**F.r*p*np.polynomial.polynomial.polyval(p, co)


def V(k, p0, kap, T, n=400, kap0=None, T1=None):
    """int_0^T |psi_k(p0 + i y)| e^{-kap y} dy  (Gauss-Laguerre-like: GL on [0, T] in the variable s = 1 - e^{-kap y})."""
    # substitution y = -log(1 - s (1 - e^{-kap T}))/kap  maps [0, 1] -> [0, T]; dy = (1 - e^{-kap T}) ds/(kap (1 - s(1 - e^{-kap T})))
    a = 1 - np.exp(-kap*T)
    x, w = np.polynomial.legendre.leggauss(n); s = (x + 1)/2; w = w/2
    y = -np.log(1 - s*a)/kap
    jac = a/(kap*(1 - s*a))
    v = float((np.abs(psi_c(k, p0 + 1j*y))*np.exp(-kap*y)*jac*w).sum())
    if kap0 is not None and T1 is not None and T1 < T:
        # rate kap on [0, T1] (= the whole integral above with T -> T1), rate kap0 < kap on [T1, T]: continuous majorant
        v = V(k, p0, kap, T1, n)
        yy = T1 + (T - T1)*s
        v += float((np.abs(psi_c(k, p0 + 1j*yy))*np.exp(-kap*T1 - kap0*(yy - T1))*w*(T - T1)).sum())
    return v


def top(k, pw, T, kapmin, Mfun, n=600):
    x, w = np.polynomial.legendre.leggauss(n); p = pw + (1 - pw)*(x + 1)/2; w = w*(1 - pw)/2
    return float(np.exp(-kapmin*T)*(np.abs(psi_c(k, p + 1j*T))*Mfun(p)*w).sum())


def kappa(beta, mfl, p, mode):
    k0 = 0.97*beta*C.L
    if mode == 'uniform':
        return k0*np.ones_like(np.asarray(p, float))
    x = np.minimum(np.asarray(p, float)/mfl, 0.999999)
    return k0/np.sqrt(1 - x*x)


def A_curv(beta, J, mfl):
    """amplitude at p = 1 of the hybrid's leading-Regge curvature waves (P_nu(cosh beta) - J0(bhat p))_+- (uniform rates Lam beta, bhat):
    1.2 x [max_{p in [0.9, 1]} |P - J0| + (beta^2/2) max|P|]  (real-axis envelope of the difference plus an allowance for the rate
    mismatch (bhat - Lam beta) y <~ beta^2/4 over the edge segment; the two waves nearly cancel at small beta)."""
    from scipy.special import j0 as _j0
    p = np.linspace(0.9, 1.0, 801); nu = C.L*p
    P = C.Pd(beta, nu)[0].real
    d = np.abs(P - _j0((2*J + 1)/mfl*p))
    return 1.2*(float(d.max()) + 0.5*beta**2*float(np.abs(P).max()))


def contour_bound(beta, mfl, pw, Mw, Mfun, mode='uniform', k=2, u=1.0, Aunif=0.0):
    """bound on |cost_c| (k = 2: 2 int psi2 G; k = 3: u int |psi3| G) given the boundary amplitude Mw at p_w and the envelope Mfun(p)
    on [p_w, 1] (used at p = 1 and on the top side).  Returns (bd, edge, top) separately.
    Aunif: amplitude at p = 1 of the part of G with only the UNIFORM rate (the hybrid's known curvature waves -(P - J0)); with mode
    'local', G = G_loc + G_unif, G_loc = exact - flat (local rate, envelope Mfun(1) + Aunif), G_unif = -(P - J0) (rate 0.97 Lam beta)."""
    w = 2.0 if k == 2 else u
    kw = float(kappa(beta, mfl, pw, mode)); k1 = float(kappa(beta, mfl, 1.0, mode)); kmin = 0.97*beta*C.L
    T = 25.0/kmin                       # AG' rectangle height: 25/(0.97 beta) in nu units (HLJ lj_scan probed y <= 30/beta)
    bd = 2*w*V(k, pw, kw, T, kap0=kmin, T1=min(T, T1F/kw))*Mw
    M1 = float(Mfun(np.array([1.0]))[0])
    if mode == 'local' and Aunif > 0:
        ed = 2*w*(V(k, 1.0, k1, T, kap0=kmin, T1=min(T, T1F/k1))*(M1 + Aunif) + V(k, 1.0, kmin, T)*Aunif)
    else:
        ed = 2*w*V(k, 1.0, k1, T, kap0=kmin, T1=min(T, T1F/k1))*(M1 + Aunif)
    tp = 2*w*top(k, pw, T, kmin, Mfun)
    return bd, ed, tp


def e_sc(beta, mfl, a0=10.0, a_bulk=4.0, a_edge=5.0, r=3.0, hel=0.0):
    """scalar-calibrated AG envelope of the beyond-hybrid correction G (per unit leading image), as a function of p:
    |G| <= [a0 + (nu beta)^2 (a_bulk + a_edge (1 - x^2)^-r)]/m^2,  nu = Lam p, x = nu/m  (fit to the exact scalar, p4w_check, x <= 0.97).
    hel: extra helicity transfer factor (1 - x^2)^hel of the hybrid envelope relative to the scalar E_2 (sm_exact / Omega^U convention:
    ++ : +2, +- : -2, off-diagonal: 0)."""
    def f(p):
        p = np.asarray(p, float); x = np.minimum(p/mfl, 0.999)
        return (a0 + (C.L*p*beta)**2*(a_bulk + a_edge*(1 - x*x)**(-r)))/(mfl*C.L)**2*(1 - x*x)**hel
    return f
