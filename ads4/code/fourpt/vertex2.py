"""
Exact frame integrand of the Lorentzian shadow representation of the doubly-spinning partial wave and its heavy limit.

Frame (CHMRS (C.2)): (x_1, x_2, x_3, x_4; x_5) = (-y, x/x^2, y, -x/x^2; e), n = (1,1,0), x, y future timelike.
For coefficient vectors c_L, c_R on the words (structures2.words), the partial wave W[c_L, c_R] has the shadow representation
   W = ((-1)^J/4) int d^3x_5 Dn  L[c_L](1,2,5) * R[c_R](3,4,5),
   L[c] = sum_w c_w D_w S3(1,2,5; Delta, J),     R[c] = sum_w c_w(Delta,J)/b_w D_w S3^S(3,4,5; 3-Delta, -1-J),
b_w = b^{(Delta_34 of the seed)} (e2e: [DVE] for scalar seeds; the D-operators commute with the x_5 integral).
For the Weyl-covariant basis v_a of t1_structures_l1 (shadow map = (1/b^{(0)}) x Weyl reflection, exact):
   R[v_a(Delta,J)] = (1/b^{(0)}) L[v_a(3-Delta,-1-J)](3,4,5; 3-Delta,-1-J).
Stripped integrand  T(z1,z2',z3,z4'; x, y) := s13^{l+1} s24^{l+1} ((-1)^J/4) L R,  with z2 = I(xhat) z2', z4 = I(xhat) z4'
(the polarisations of operators 2, 4 are referred to xhat by the inversion tensor).
Moments:  Pihat^q[s13^{l+1} s24^{l+1} W] = 2 int d^3x d^3y (|x||y|)^{k-4} K^q(xhat, yhat) T(x, y)   (e2e eq. (*)).
"""
import os, sys, functools
import mpmath as mp, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
from formulas import log_b                       # noqa: E402
from structures2 import rel_structures, logMref, words
from embed2 import Delta as sD, J as sJ

def mdot(u, v): return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def msub(u, v): return [u[i] - v[i] for i in range(3)]
def madd(u, v): return [u[i] + v[i] for i in range(3)]
def msc(a, u): return [a*u[i] for i in range(3)]

@functools.lru_cache(None)
def struct_funs(l):
    rs, S, v5 = rel_structures(l)
    names = ['s12', 's15', 's25', 'z12', 'z15', 'z21', 'z25', 'z51', 'z52', 'z1z2', 'z1z5', 'z2z5']
    args = [S[n] for n in names] + [sD, sJ]
    fq = {w: sp.lambdify(args, Q, modules='mpmath') for (w, Q) in rs}
    lm = logMref(l, S, v5)
    fl = sp.lambdify(args, lm.subs(sp.log(v5), sp.log(sp.Abs(v5))), modules='mpmath')
    return fq, fl

def structure(l, coeffs, D, Jv, xa, xb, x5, za, zb, n):
    """sum_w coeffs[w] W_w at points (a, b, 5): returns (log of positive prefactor M_ref, polynomial part)."""
    fq, fl = struct_funs(l)
    sab = -mdot(msub(xa, xb), msub(xa, xb)); sa5 = -mdot(msub(xa, x5), msub(xa, x5)); sb5 = -mdot(msub(xb, x5), msub(xb, x5))
    a = [sab, sa5, sb5, mdot(za, msub(xb, xa)), mdot(za, msub(x5, xa)), mdot(zb, msub(xa, xb)), mdot(zb, msub(x5, xb)),
         mdot(n, msub(xa, x5)), mdot(n, msub(xb, x5)), mdot(za, zb), mdot(za, n), mdot(zb, n), D, Jv]
    return fl(*a), sum(coeffs[w]*fq[w](*a) for w in coeffs)

def inv(x, v):
    """inversion tensor I(x) v = v - 2 x (x.v)/x^2"""
    x2 = mdot(x, x)
    return msub(v, msc(2*mdot(x, v)/x2, x))

def integrand(l, cL, cR_refl, D, Jv, x, y, z1, z2p, z3, z4p, logb=None):
    """T for left coefficients cL (at (Delta,J)) and right coefficients cR_refl = c_R evaluated at (3-Delta,-1-J) times
    b^{(0)}/b_w (i.e. the right structure is (1/b^{(0)}) L[cR_refl](3,4,5; 3-Delta, -1-J)).  Returns (logscale, value):
    T = exp(logscale)*value  (the sign (-1)^J included in value)."""
    D = mp.mpf(D)
    e = [mp.mpf(1), mp.mpf(0), mp.mpf(0)]; n = [mp.mpf(1), mp.mpf(1), mp.mpf(0)]
    x2n = mdot(x, x)
    x1 = [-v for v in y]; x2 = [v/x2n for v in x]; x3 = list(y); x4 = [-v/x2n for v in x]
    z2 = inv(x, z2p); z4 = inv(x, z4p)
    lL, pL = structure(l, cL, D, Jv, x1, x2, e, z1, z2, n)
    lR, pR = structure(l, cR_refl, 3 - D, -1 - Jv, x3, x4, e, z3, z4, n)
    s13 = -mdot(msub(x1, x3), msub(x1, x3)); s24 = -mdot(msub(x2, x4), msub(x2, x4))
    if logb is None:
        logb = log_b(D, Jv, 0)
    logscale = lL + lR + (l+1)*(mp.log(s13) + mp.log(s24)) - mp.log(4) - logb
    return logscale, (-1)**Jv*pL*pR

# ---------------------------------------------------------------- frames on H^2 (rule: e_par at a points to b; e_perp.v = det[v, a, e_par])
def det3(a, b, c):
    return (a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0]) + a[2]*(b[0]*c[1] - b[1]*c[0]))

def frame(a, b):
    """unit tangent at a (in H^2) pointing to b, and e_perp with e_perp.v = det[v, a, e_par] (contravariant components)."""
    ch = -mdot(a, b)
    sh = mp.sqrt(ch*ch - 1)
    epar = msc(1/sh, msub(b, msc(ch, a)))
    # e_perp: covariant components (det with unit vectors), raise index with eta = diag(-1,1,1)
    cov = [det3([1 if i == k else 0 for i in range(3)], a, epar) for k in range(3)]
    eperp = [-cov[0], cov[1], cov[2]]
    return epar, eperp

def eps(epar, eperp, h):
    s = mp.sqrt(2)
    return [(epar[i] + 1j*h*eperp[i])/s for i in range(3)]

def helicity_components(l, cL, cR_refl, D, Jv, X, Y, m_scale):
    """heavy-limit helicity coefficients at fixed (X, Y) = (m x, m y):  T / [((-1)^J/(4 b)) 2 4^{2(l+1)} |x|^{2(l+1)} |y|^{2(l+1)}
    e^{2P.x + 2Pbar.y}],  P = (A_+, A_-, 0), Pbar = (A_+, -A_-, 0)  evaluated on z1 = conj eps_{h1}(pbarhat) etc. (l = 1).
    Returns dict (h1,h2,h3,h4) -> value, plus the longitudinal leakage."""
    assert l == 1
    D = mp.mpf(D)
    Ap, Am = D - mp.mpf(3)/2, Jv + mp.mpf(1)/2
    m = mp.sqrt(Ap*Ap - Am*Am)
    x = [v/m_scale for v in X]; y = [v/m_scale for v in Y]
    P = [Ap, Am, mp.mpf(0)]; Pb = [Ap, -Am, mp.mpf(0)]
    ph = msc(1/m, P); pbh = msc(1/m, Pb)
    fL = frame(pbh, ph); fR = frame(ph, pbh)
    ax = mp.sqrt(-mdot(x, x)); ay = mp.sqrt(-mdot(y, y))
    lpw = (mp.log(2) + 2*(l+1)*mp.log(4) + 2*(l+1)*(mp.log(ax) + mp.log(ay)) + 2*mdot(P, x) + 2*mdot(Pb, y) - mp.log(4) - log_b(D, Jv, 0))
    res = {}
    for h1 in (1, -1):
        for h2 in (1, -1):
            for h3 in (1, -1):
                for h4 in (1, -1):
                    z1 = eps(*fL, -h1); z3 = eps(*fL, -h3); z2p = eps(*fR, -h2); z4p = eps(*fR, -h4)
                    ls, v = integrand(l, cL, cR_refl, D, Jv, x, y, z1, z2p, z3, z4p)
                    res[(h1, h2, h3, h4)] = v*mp.exp(ls - lpw)*(-1)**Jv
    # leakage: z1 along pbarhat
    ls, v = integrand(l, cL, cR_refl, D, Jv, x, y, pbh, eps(*fR, 1), eps(*fL, 1), eps(*fR, -1))
    leak = v*mp.exp(ls - lpw)*(-1)**Jv
    return res, leak, m
