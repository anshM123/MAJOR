"""
Hyperbolic (R x H^2) shear response engine: Einstein + theta(r) I and theta(r) J on the AdS4 hyperbolic black hole.

Background (L = R_H2 = 1):  ds^2 = -f dt^2 + dr^2/f + r^2 gamma,  f = r^2 - 1 - mu/r,  mu = rh^3 - rh,
H^2 in coordinates u = tanh X:  gamma = du^2/w^2 + dY^2/w,  w = 1 - u^2   (Gauss curvature -1).
TT tensor on H^2: H = h (dX^2 - cosh^2X dY^2), h = sech^2 X = w  (divergence-free and traceless; L^2 in X).
Perturbation (det-preserving, as in ../../scope_thermal/allorders/static_engine.py):
   g_uu = r^2 e^{eps phi(r) w}/w^2,   g_YY = r^2 e^{-eps phi(r) w}/w.
The O(eps^2) Lagrangian density is integrated over u in (-1,1) with sqrt(-g) = r^2 w^{-3/2} (Y is a spectator), then reduced to
sum_k alpha_k(r) (phi^(k))^2 and eta is read off with the static membrane flux of static_engine (horizon at r = rh).
Horizon Weyl scale: x = mu/rh^3 = 1 - 1/rh^2  (x = 0: AdS-Rindler = vacuum at T0 = 1/2pi;  x -> 1: planar brane).
"""
import sys, itertools
import sympy as sp
sys.path.insert(0, '/home/user/MAJOR/ads4/code/scope_thermal/allorders')
import static_engine as se
from static_engine import r, eps, P, phi, eps_series, reduce_quadratic

u = se.X[2]                       # use the engine's x-coordinate as u
rh = sp.symbols('r_h', positive=True)
w = 1 - u**2


def background_f():
    mu = rh**3 - rh
    return sp.factor(r**2 - 1 - mu/r)


def hyp_metric(f, pert=True):
    e = sp.exp(eps*phi*w) if pert else 1
    return sp.diag(-f, 1/f, r**2*e/w**2, r**2/(e*w))


def integrate_u(expr):
    """Int_{-1}^{1} expr * w^{-3/2} du for expr rational in u (even part), via w-powers and Beta functions."""
    expr = sp.together(sp.expand(expr))
    num, den = sp.fraction(expr)
    # den must be (possibly) r-stuff times powers of w
    den_f = sp.factor(den)
    k = 0
    rest = den_f
    while True:
        q, rem = sp.div(sp.Poly(rest, u), sp.Poly(w, u))
        if rem.is_zero and q.degree() >= 0 and sp.Poly(rest, u).degree() >= 2:
            rest = q.as_expr(); k += 1
        else:
            break
    if sp.Poly(rest, u).degree() > 0:
        raise ValueError('unexpected u-dependence in denominator: %s' % rest)
    pn = sp.Poly(sp.expand(num), u)
    tot = 0
    for (deg,), c in pn.terms():
        if deg % 2:
            continue
        # u^deg = (1-w)^(deg/2)
        m = deg//2
        for j in range(m + 1):
            cj = sp.binomial(m, j)*(-1)**j      # coefficient of w^j
            pw = j - k - sp.Rational(3, 2)
            if pw <= -1:
                tot += sp.Symbol('DIVERGENT_%d' % j)*c*cj
            else:
                tot += c*cj*sp.sqrt(sp.pi)*sp.gamma(pw + 1)/sp.gamma(pw + sp.Rational(3, 2))
    return sp.simplify(tot/rest)


def reduced_lagrangian(density2):
    """density2: O(eps^2) coefficient (jets p0..p4, r, u) of the scalar density WITHOUT sqrt(-g); returns sum c_ab p_a p_b
    after the u-integration (sqrt(-g) = r^2 w^{-3/2} included)."""
    e = sp.expand(density2*r**2)
    out = 0
    for a in range(5):
        for b in range(a, 5):
            if a == b:
                c = e.coeff(P[a], 2)
                for o in range(5):
                    if o != a:
                        c = c.subs(P[o], 0)
            else:
                c = e.coeff(P[a], 1).coeff(P[b], 1)
            if c != 0:
                out += integrate_u(c)*P[a]*P[b]
    return out


def flux_at_horizon(alpha, rsp):
    chi_der = {1: -rsp}
    for k in range(2, 5):
        chi_der[k] = sp.diff(chi_der[k - 1], r)
    Pi = 0
    for k, a in alpha.items():
        if k == 0:
            continue
        term = 2*a*chi_der[k]
        for _ in range(k - 1):
            term = sp.diff(term, r)
        Pi += (-1)**(k - 1)*term
    return sp.limit(sp.simplify(-Pi), r, rh)
