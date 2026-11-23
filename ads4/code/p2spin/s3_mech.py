"""
p2spin/s3_mech.py (P2_SPINNING.md §3.1): CHMRS's p-derivative mechanism (2106.10274 §3.3) as an exact symbolic map
    polynomial insertion f(X, Y)  ->  ring-class coefficients at the state's beta,
and its validation on the scalar seed against CHMRS's ancillary series (regge_moments_order8.m, 'series Pi', order 1/m2).

Mechanism.  With the leading symmetric-frame integrand  |X|^{a-3}|Y|^{a-3} P_nu(eta_XY) exp(2 P.X + 2 Pb.Y)  and the
Fourier-Laplace identity (CHMRS eq. 'Fourier identity') applied twice, for independent P, Pb:
    F(P, Pb) = |P|^-a |Pb|^-a  P_nu(-P.Pb/(|P||Pb|))    (up to nu-only gamma factors),
an insertion X^al Y^be becomes (1/2 d_P)^al (1/2 d_Pb)^be F.  Higher derivatives of P_nu are reduced with Legendre's equation
(eta^2 - 1) P'' + 2 eta P' + (nu^2 + 1/4) P = 0 (d = 3).  The result, divided by F at the state, is c0(eta, nu^2) P + c1(eta, nu^2) P'.

Seed: CHMRS's T(x, y) (eq. 'moment appendix'), scalar Delta_phi = 3, d = 3:
    T = (x^2 y^2)^{-tD} (1 + 2x.y + x^2y^2)^{d-1-2Dphi} / ((e-x)^2 (e-y)^2) * R^{(Delta+J)/2} * [-n.V_e]^J / [n.V_{-e}]^{J+d-2},
expanded at x, y = eps X, eps Y with Delta + J - 1 = q/eps, Delta - J - 2 = 1/(q eps).
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp

EPS = sp.Symbol('eps')
nu2, eta = sp.symbols('nu2 eta')                    # nu^2, eta = cosh beta
Pf, dPf = sp.symbols('Pf dPf')                      # P_nu(eta), P_nu'(eta)
X = sp.symbols('X0:3'); Y = sp.symbols('Y0:3')      # Cartesian (t, z, y), metric (-,+,+)
MET = sp.diag(-1, 1, 1)


def mdot(a, b):
    return sum(MET[i, i]*a[i]*b[i] for i in range(3))


# ----------------------------------------------------------------------------------------------- seed expansion
def seed_series(qv, Dphi=3, d=3, order=2):
    """log T expanded: returns (Pvec, Pbvec, [s1, s2]) with T = const*exp(2P.X + 2Pb.Y)*(1 + eps s1 + eps^2 s2)*(X^2Y^2)^{-tD}."""
    e = (1, 0, 0); n = (1, 1, 0)
    x = [EPS*c for c in X]; y = [EPS*c for c in Y]
    yb = [-y[i] - 2*e[i]*mdot(y, e) for i in range(3)]
    D = ((qv + 1/qv)/EPS + 3)/2; J = ((qv - 1/qv)/EPS - 1)/2
    def V(sg):
        ee = [sg*c for c in e]
        a = (1 - 2*mdot(ee, x))*(1 - 2*mdot(ee, y)) - mdot(x, x)*mdot(y, y)
        ex2 = mdot([ee[i] + x[i] for i in range(3)], [ee[i] + x[i] for i in range(3)])
        ey2 = mdot([ee[i] + y[i] for i in range(3)], [ee[i] + y[i] for i in range(3)])
        ybs = [-y[i] - 2*ee[i]*mdot(y, ee) for i in range(3)] if sg == 1 else [-y[i] - 2*e[i]*mdot(y, e) for i in range(3)]
        return [a*ee[i] + ey2*x[i] + ex2*ybs[i] for i in range(3)]
    sq = lambda u: mdot(u, u)
    em = lambda u, s: [e[i] + s*u[i] for i in range(3)]
    terms = [((d - 1 - 2*Dphi), 1 + 2*mdot(x, y) + sq(x)*sq(y)),
             (-1, sq(em(x, -1))), (-1, sq(em(y, -1))),
             ((D + J)/2, sq(em(x, -1))*sq(em(y, -1))/(sq(em(x, 1))*sq(em(y, 1)))),
             (J, -mdot(n, V(1))), (-(J + d - 2), mdot(n, V(-1)))]
    L = 0
    for pw, arg in terms:
        a0 = sp.expand(arg.subs(EPS, 0))
        L += pw*sp.series(sp.log(sp.expand(arg/a0)), EPS, 0, order + 2).removeO()   # constant log a0 dropped
    L = sp.expand(sp.series(sp.expand(L), EPS, 0, order + 1).removeO())
    L0 = L.coeff(EPS, 0); L1 = L.coeff(EPS, 1); L2 = L.coeff(EPS, 2)
    P = [sp.Rational(1, 2)*sp.diff(L0, X[i])*MET[i, i] for i in range(3)]
    Pb = [sp.Rational(1, 2)*sp.diff(L0, Y[i])*MET[i, i] for i in range(3)]
    s1 = sp.expand(L1); s2 = sp.expand(L2 + L1**2/2)
    return P, Pb, [s1, s2]


# ----------------------------------------------------------------------------------------------- the mechanism
_cache = {}


def deriv_table(a, maxdeg=6):
    """Symbolic F and its P-, Pb-derivatives with generic P, Pb (Cartesian symbols); F = |P|^-a |Pb|^-a G(eta(P, Pb))."""
    Ps = sp.symbols('p0:3'); Qs = sp.symbols('pb0:3')
    G = sp.Function('G')
    nP = sp.sqrt(-mdot(Ps, Ps)); nQ = sp.sqrt(-mdot(Qs, Qs))
    et = -mdot(Ps, Qs)/(nP*nQ)
    F = nP**(-a)*nQ**(-a)*G(et)
    return Ps, Qs, G, F


def apply_insertion(poly, P0, Pb0, a):
    """I[poly]/F at (P0, Pb0) as c0*Pf + c1*dPf, exact (poly: polynomial in X, Y)."""
    Ps, Qs, G, F = deriv_table(a)
    poly = sp.Poly(sp.expand(poly), *X, *Y)
    sub = dict(zip(Ps, P0)); sub.update(dict(zip(Qs, Pb0)))
    F0 = F.subs(sub)
    tot = 0
    for mon, c in poly.terms():
        expr = F
        for i in range(3):
            for _ in range(mon[i]):
                expr = sp.diff(expr, Ps[i])*sp.Rational(1, 2)*MET[i, i]      # X^i insertion = (1/2) d/dP_i (upper index)
        for i in range(3):
            for _ in range(mon[3 + i]):
                expr = sp.diff(expr, Qs[i])*sp.Rational(1, 2)*MET[i, i]
        tot += c*expr.subs(sub)
    # replace G and its derivatives at eta0
    et0 = sp.nsimplify(sp.simplify((-mdot(P0, Pb0))/sp.sqrt(mdot(P0, P0)*mdot(Pb0, Pb0))))
    tot = sp.simplify(tot/F0)
    w = sp.Symbol('w')
    reps = {}
    for k in range(6, -1, -1):
        reps[k] = sp.Symbol(f'G{k}')
    tot = tot.replace(lambda e: isinstance(e, sp.Subs), lambda e: reps[sum(1 for _ in e.expr.variables) if hasattr(e.expr, 'variables') else 0])
    tot = tot.subs({sp.Derivative(G(w), (w, k)).subs(w, et0): reps[k] for k in range(1, 7)})
    tot = tot.subs(G(et0), reps[0])
    return sp.expand(tot), et0, reps


def reduce_legendre(expr, reps, et0):
    """G_k -> c0 G0 + c1 G1 via (e^2-1)G'' + 2e G' + (nu2 + 1/4) G = 0, then G0 -> Pf/Pf (normalised), G1 -> dPf."""
    e = et0
    # G'' = -(2 e G' + (nu2 + 1/4) G)/(e^2 - 1); higher by differentiating the ODE: (e^2-1)G^{(k+2)} + 2(k+1) e G^{(k+1)} + (k(k+1) + nu2 + 1/4) G^{(k)} = 0
    sol = {0: (1, 0), 1: (0, 1)}
    for k in range(0, 5):
        c0 = -(2*(k + 1)*e*sol[k + 1][0] + (k*(k + 1) + nu2 + sp.Rational(1, 4))*sol[k][0])/(e**2 - 1)
        c1 = -(2*(k + 1)*e*sol[k + 1][1] + (k*(k + 1) + nu2 + sp.Rational(1, 4))*sol[k][1])/(e**2 - 1)
        sol[k + 2] = (sp.simplify(c0), sp.simplify(c1))
    out0 = 0; out1 = 0
    for k, s in reps.items():
        cf = sp.expand(expr).coeff(s)
        out0 += cf*sol[k][0]; out1 += cf*sol[k][1]
    return sp.factor(sp.simplify(out0)), sp.factor(sp.simplify(out1))


def chmrs_seriesPi(d=3, Dphi=3, alpha=None):
    """CHMRS ancillary 'series Pi' at order 1/m2, per unit leading P (divided by 2): (coef of P, coef of P')."""
    al = alpha
    cP = eta*((8 - 8*d - 6*d**2 + 10*d**3 + 4*al - 12*d*al + 3*d**2*al - 12*d*al**2 - 4*al**3 - 8*Dphi + 24*d*Dphi
               - 54*d**2*Dphi + 24*al**2*Dphi + 96*d*Dphi**2 - 64*Dphi**3)/12 + (2 + 2*d + 3*al - 6*Dphi)*nu2/3)
    cD = ((-12 + 36*d - 13*d**2 + 24*al + 12*d*al + 12*al**2 - 48*Dphi + 24*d*Dphi - 48*al*Dphi)/12 - nu2/3
          + eta**2*((-20 - 12*d + 9*d**2 - 24*al - 12*d*al - 12*al**2 + 48*Dphi - 24*d*Dphi + 48*al*Dphi)/12 + nu2/3))
    return sp.expand(cP/2), sp.expand(cD/2)
