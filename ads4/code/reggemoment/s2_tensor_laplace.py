"""
Stage 2(b): Fourier-Laplace transform of tensor ("spinning") Regge-moment kernels in d = 3.

Scalar input (CHMRS (3.18)/(B.2), d = 3; tested numerically in part (0) below):
   int_{y>0} d^3y |y|^{b-3} P_{-1/2+i nu}(-yhat.xhat) e^{y.p} = C_b gamma_b(nu) Psi_b(p),
   Psi_b(p) = |p|^{-b} P_{-1/2+i nu}(-phat.xhat),  C_b = 2^{b-1} sqrt(pi),  gamma_b(nu) = |Gamma((b-1/2+i nu)/2)|^2.

Tensor kernels at the point yhat of H^2 (operators 1 and 3 of the Regge pair; z1, z3 their polarisations), acting
on Phi_a(y) = |y|^{a-3} P(-yhat.xhat)   [structures of Costa-Hansen-Penedones (2.41), plus a parity-odd one]:
   D1 = z1.z3 + (z1.yhat)(z3.yhat)                     (transverse metric at yhat)
   D2 = -(z1.yhat)(z3.yhat)
   D3 = (z1.yhat)(z3.nabla) + (z3.yhat)(z1.nabla)
   D4 = -(z1.nabla)(z3.nabla) + D1 nabla^2/2
   D5 = eps(z1, z3, yhat)                               (parity odd)
Each D_s Phi_a is rewritten as (polynomial in y) x (flat derivatives) of Phi_b's [checked in coordinates, part (1)],
so that its Laplace transform is (polynomial in d/dp) x (powers of -p) acting on C_b gamma_b Psi_b.  The result is
then contracted with the transverse helicity vectors at phat,
   eps_pm = (n_par +- i e_perp)/sqrt2,   n_par = unit tangent of the geodesic xhat -> phat at phat (pointing away from xhat),
and expressed through P(cosh xi), P^{(2)} := sinh^2 xi P''(cosh xi)  (xi = dist(xhat, phat)):
   helicity-preserving  (z1, z3) = (eps_+, eps_-):   L[D_s Phi_a] = C_a gamma_a |p|^{-a} * M0_s(nu, a) * P
   helicity-flip        (z1, z3) = (eps_+, eps_+):   L[D_s Phi_a] = C_a gamma_a |p|^{-a} * M2_s(nu, a) * P^{(2)}
Part (2) derives M0_s, M2_s symbolically; part (3) checks them by direct numerical integration;
part (4): stress-tensor pair, non-derivative structures.
"""
import sympy as sp, mpmath as mp, itertools

nu, a = sp.symbols('nu a', positive=True)
lam = nu**2 + sp.Rational(1, 4)
ETA = sp.diag(-1, 1, 1)
def dot(u, v):
    return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def eps3(u, v, w):
    return sp.Matrix([list(u), list(v), list(w)]).det()

Pf = sp.Function('P')
c = sp.Symbol('c')
def reduce_P(expr, cc):
    """Use Legendre's equation (c^2-1)P'' + 2cP' - lam... i.e. P'' = -(2 c P' + lam P)/(c^2-1) to reduce to P, P'."""
    P0, P1 = sp.symbols('P0 P1')
    d2 = -(2*c*P1 + lam*P0)/(c**2 - 1)
    d3 = sp.diff(d2, c) + sp.diff(d2, P0)*P1 + sp.diff(d2, P1)*d2
    d4 = sp.diff(d3, c) + sp.diff(d3, P0)*P1 + sp.diff(d3, P1)*d2
    expr = expr.subs(sp.Derivative(Pf(c), (c, 4)), d4).subs(sp.Derivative(Pf(c), (c, 3)), d3)
    expr = expr.subs(sp.Derivative(Pf(c), (c, 2)), d2).subs(sp.Derivative(Pf(c), c), P1).subs(Pf(c), P0)
    return expr.subs(c, cc), P0, P1

# ------------------------------------------------------------------ generic machinery in Cartesian coordinates
v = sp.symbols('v0:3', real=True)             # generic point of the forward cone (used both as y and as p)
norm = sp.Symbol('R', positive=True)             # |v| = sqrt(v0^2 - v1^2 - v2^2), kept as a symbol; chain rule by hand
NORM_EXPL = sp.sqrt(v[0]**2 - v[1]**2 - v[2]**2)
cv = v[0]/norm                                 # = -vhat.xhat for xhat = (1,0,0)
_dR = [v[0]/norm, -v[1]/norm, -v[2]/norm]
_dc = [1/norm - v[0]**2/norm**3, v[0]*v[1]/norm**3, v[0]*v[2]/norm**3]
def Phi(b):                                    # |v|^{b-3} P(cv)
    return norm**(b-3)*Pf(c)
def Psi(b):                                    # |v|^{-b} P(cv)
    return norm**(-b)*Pf(c)
def D(expr, mu):                               # d/dv^mu, chain rule through c and R
    return sp.diff(expr, v[mu]) + sp.diff(expr, c)*_dc[mu] + sp.diff(expr, norm)*_dR[mu]
def zD(z, expr):                               # z^mu d/dv^mu
    return sum(z[m]*D(expr, m) for m in range(3))
def zv(z):
    return dot(z, v)

def structures_y(z1, z3, aa):
    """D_s Phi_a written with flat derivatives (the representation used for the Laplace transform)."""
    out = {}
    out[1] = dot(z1, z3)*Phi(aa) + zv(z1)*zv(z3)*Phi(aa-2)
    out[2] = -zv(z1)*zv(z3)*Phi(aa-2)
    out[3] = zv(z1)*zD(z3, Phi(aa)) + zv(z3)*zD(z1, Phi(aa)) + 2*(aa-3)*zv(z1)*zv(z3)*Phi(aa-2)
    hess = (zD(z1, zD(z3, Phi(aa+2))) + (aa-2)*(zv(z1)*zD(z3, Phi(aa)) + zv(z3)*zD(z1, Phi(aa)))
            + (aa-1)*dot(z1, z3)*Phi(aa) + (aa-3)**2*zv(z1)*zv(z3)*Phi(aa-2))
    out[4] = -hess - lam/2*out[1]
    out[5] = eps3(z1, z3, v)*Phi(aa-1)
    return out

def structures_direct(z1, z3, aa):
    """The same structures from their definition: tangential projections and the covariant Hessian on H^2,
    evaluated for the degree-0 extension of P (independent code path for the check in part (1))."""
    vh = [x/norm for x in v]
    P0 = Pf(c)
    def tang(z):                                 # projection on the tangent space at vhat
        return [z[m] + vh[m]*dot(vh, z) for m in range(3)]
    def grad(f):                                 # flat gradient with upper index
        return [sum(ETA[m, n]*D(f, n) for n in range(3)) for m in range(3)]
    z1t, z3t = tang(z1), tang(z3)
    g = grad(P0)
    nab = lambda zt: norm*sum(zt[m]*D(P0, m) for m in range(3))               # z.nabla P  (unit hyperboloid)
    hess = norm**2*sum(z1t[m]*z3t[n]*D(D(P0, n), m) for m in range(3) for n in range(3))   # tangential Hessian of a degree-0 function
    D1 = dot(z1, z3) + dot(z1, vh)*dot(z3, vh)
    out = {1: D1*P0, 2: -dot(z1, vh)*dot(z3, vh)*P0, 3: dot(z1, vh)*nab(z3t) + dot(z3, vh)*nab(z1t),
           4: -hess - lam/2*D1*P0, 5: eps3(z1, z3, vh)*P0}
    return {s: norm**(aa-3)*e for s, e in out.items()}

def laplace_structures(z1, z3, aa):
    """Laplace transforms: y^mu -> d/dp_mu, d/dy^mu -> -p_mu, Phi_b -> C_b gamma_b Psi_b.  Returned divided by
    C_a gamma_a, using C_{b}/C_a = 2^{b-a} and gamma_{a-2}/gamma_a = 4/((a-5/2)^2+nu^2), gamma_{a+2}/gamma_a = ((a-1/2)^2+nu^2)/4,
    gamma_{a-1} kept as the symbol g1 := gamma_{a-1}/gamma_a."""
    g1 = sp.Symbol('g1')
    rel = {aa: 1, aa-2: sp.Rational(1, 4)*4/((aa-sp.Rational(5, 2))**2 + nu**2),
           aa+2: 4*((aa-sp.Rational(1, 2))**2 + nu**2)/4, aa-1: sp.Rational(1, 2)*g1}
    T = lambda b: rel[b]*Psi(b)
    zp = lambda z: dot(z, v)                     # z.p
    out = {}
    out[1] = dot(z1, z3)*T(aa) + zD(z1, zD(z3, T(aa-2)))
    out[2] = -zD(z1, zD(z3, T(aa-2)))
    out[3] = zD(z1, -zp(z3)*T(aa)) + zD(z3, -zp(z1)*T(aa)) + 2*(aa-3)*zD(z1, zD(z3, T(aa-2)))
    hess = (zp(z1)*zp(z3)*T(aa+2) + (aa-2)*(zD(z1, -zp(z3)*T(aa)) + zD(z3, -zp(z1)*T(aa)))
            + (aa-1)*dot(z1, z3)*T(aa) + (aa-3)**2*zD(z1, zD(z3, T(aa-2))))
    out[4] = -hess - lam/2*out[1]
    out[5] = sum(sp.diff(eps3(z1, z3, v), v[m])*sum(ETA[m, n]*D(T(aa-1), n) for n in range(3)) for m in range(3))
    return out

xi = sp.Symbol('xi', positive=True)
tt = sp.Symbol('t', positive=True)                               # t = tanh(xi/2): cosh xi = (1+t^2)/(1-t^2), sinh xi = 2t/(1-t^2)
CH, SH = (1 + tt**2)/(1 - tt**2), 2*tt/(1 - tt**2)
POINT = {v[0]: CH, v[1]: SH, v[2]: 0}                            # |p| = 1 (the power |p|^{-a} is restored by homogeneity)
N_PAR = [SH, CH, 0]
E_PERP = [0, 0, 1]
EP = [(N_PAR[m] + sp.I*E_PERP[m])/sp.sqrt(2) for m in range(3)]
EM = [(N_PAR[m] - sp.I*E_PERP[m])/sp.sqrt(2) for m in range(3)]

def at_point(expr):
    """evaluate at p = (cosh xi, sinh xi, 0), |p| = 1: everything becomes rational in t = tanh(xi/2)."""
    e = expr.subs(norm, 1).subs(POINT)
    e, P0, P1 = reduce_P(e, CH)
    return sp.cancel(sp.together(e)), P0, P1

def derive_M(verbose=True):
    res = {}
    for name, (z1, z3) in {'0(+-)': (EP, EM), '0(-+)': (EM, EP), '2(++)': (EP, EP)}.items():
        L = laplace_structures(z1, z3, a)
        for s, e in L.items():
            e, P0, P1 = at_point(e)
            if name.startswith('0'):
                coeffP1 = sp.cancel(sp.diff(e, P1)); M = sp.factor(sp.cancel(sp.diff(e, P0)))
                assert coeffP1 == 0, (name, s, coeffP1)
            else:
                # P^{(2)} = sinh^2 xi P'' = -(2 cosh xi P' + lam P)
                M = sp.factor(sp.cancel(sp.diff(e, P1)/(-2*CH)))
                assert sp.cancel(e - M*(-(2*CH*P1 + lam*P0))) == 0, (name, s)
            assert tt not in M.free_symbols, (name, s, M)
            res[(name, s)] = M
            if verbose:
                print(f"   M[{name}]_{s} = {M}")
    return res

if __name__ == "__main__":
    mp.mp.dps = 16
    print("(0) scalar transform (B.2) in d = 3, direct 2D numerical integration, p boosted relative to xhat")
    def scalar_lhs(b, nuv, xi0):
        # y = s (cosh t, sinh t cos f, sinh t sin f); int ds s^{b-1} e^{-s w} = Gamma(b) w^{-b}, w = -yhat.p
        Pn = lambda x: mp.re(mp.hyp2f1(mp.mpf(1)/2 - 1j*nuv, mp.mpf(1)/2 + 1j*nuv, 1, (1-x)/2))
        f = lambda t, ph: mp.sinh(t)*Pn(mp.cosh(t))*(mp.cosh(t)*mp.cosh(xi0) - mp.sinh(t)*mp.sinh(xi0)*mp.cos(ph))**(-b)
        return mp.gamma(b)*2*mp.quad(f, [0, 1, 3, 8, 30], [0, mp.pi])
    for (b, nuv, xi0) in ((mp.mpf('3.5'), mp.mpf('0.8'), mp.mpf('0.6')), (mp.mpf('5'), mp.mpf('2.0'), mp.mpf('1.3'))):
        Pn = mp.re(mp.hyp2f1(mp.mpf(1)/2 - 1j*nuv, mp.mpf(1)/2 + 1j*nuv, 1, (1-mp.cosh(xi0))/2))
        rhs = 2**(b-1)*mp.sqrt(mp.pi)*abs(mp.gamma((b - mp.mpf(1)/2 + 1j*nuv)/2))**2*Pn
        print(f"   b={b} nu={nuv} xi={xi0}:  lhs {mp.nstr(scalar_lhs(b, nuv, xi0), 12)}   rhs {mp.nstr(rhs, 12)}")

    print("(1) flat-derivative representation of D_s Phi_a vs definition (random numerical point, generic z1, z3)")
    import random
    random.seed(3)
    z1 = [sp.Rational(random.randint(-5, 5), 3) for _ in range(3)]; z3 = [sp.Rational(random.randint(-5, 5), 2) for _ in range(3)]
    A = structures_y(z1, z3, a); B = structures_direct(z1, z3, a)
    pt = {v[0]: sp.Rational(17, 5), v[1]: sp.Rational(-3, 4), v[2]: sp.Rational(11, 7)}
    for s in A:
        Rv = NORM_EXPL.subs(pt)
        dA, P0, P1 = reduce_P(A[s].subs(norm, Rv).subs(pt), (v[0]/Rv).subs(pt)); dB, _, _ = reduce_P(B[s].subs(norm, Rv).subs(pt), (v[0]/Rv).subs(pt))
        print(f"   D{s}: difference = {sp.simplify(dA - dB)}")

    print("(2) symbolic transforms contracted with helicity vectors at phat:  L[D_s Phi_a] = C_a gamma_a |p|^{-a} M P  resp. M P^(2);  g1 := gamma_{a-1}/gamma_a")
    M = derive_M()

    print("(3) numerical check of (2) by direct integration of the tensor kernel (a = 4.5, nu = 1.1, xi = 0.9); Gauss-Legendre, float64")
    import numpy as np
    av, nuv, xiv = 4.5, 1.1, 0.9
    _A, _B = mp.mpf(1)/2 - 1j*mp.mpf(nuv), mp.mpf(1)/2 + 1j*mp.mpf(nuv)
    def Pd(x, n):                    # d^n/dx^n 2F1(A,B;1;(1-x)/2) = (-1/2)^n (A)_n (B)_n / n! 2F1(A+n,B+n;1+n;(1-x)/2)
        return float(mp.re((-mp.mpf(1)/2)**n*mp.rf(_A, n)*mp.rf(_B, n)/mp.factorial(n)*mp.hyp2f1(_A+n, _B+n, 1+n, (1-mp.mpf(x))/2)))
    gam = lambda b: float(abs(mp.gamma((b - mp.mpf(1)/2 + 1j*mp.mpf(nuv))/2))**2)
    def glnodes(a0, b0, n):
        x, w = np.polynomial.legendre.leggauss(n)
        return 0.5*(b0-a0)*x + 0.5*(b0+a0), 0.5*(b0-a0)*w
    tn, tw = map(np.concatenate, zip(*[glnodes(a0, b0, 60) for a0, b0 in ((0, 1.5), (1.5, 4), (4, 8), (8, 13))]))
    fn, fw = glnodes(0, 2*np.pi, 160)
    P0a = np.array([Pd(np.cosh(t), 0) for t in tn]); P1a = np.array([Pd(np.cosh(t), 1) for t in tn]); P2a = np.array([Pd(np.cosh(t), 2) for t in tn])
    T, F = np.meshgrid(tn, fn, indexing='ij'); Wt = np.outer(tw, fw)
    yh = [np.cosh(T), np.sinh(T)*np.cos(F), np.sinh(T)*np.sin(F)]
    pvec = [np.cosh(xiv), np.sinh(xiv), 0.0]
    wdot = yh[0]*pvec[0] - yh[1]*pvec[1] - yh[2]*pvec[2]
    P0s, P1s, P2s = sp.symbols('P0s P1s P2s')
    def lamb(expr):
        # evaluated on the unit hyperboloid: |y| = 1
        e = expr.subs(sp.Derivative(Pf(c), (c, 2)), P2s).subs(sp.Derivative(Pf(c), c), P1s).subs(Pf(c), P0s).subs(c, cv).subs(norm, 1)
        return sp.lambdify((v[0], v[1], v[2], P0s, P1s, P2s), e.subs({a: sp.Rational(9, 2), nu: sp.Rational(11, 10)}), 'numpy')
    sub = {tt: sp.tanh(sp.Rational(9, 20))}
    from math import gamma as _G
    for name, (z1, z3) in {'0(+-)': (EP, EM), '2(++)': (EP, EP)}.items():
        z1n = [complex(sp.N(x.subs(sub), 20)) for x in z1]; z3n = [complex(sp.N(x.subs(sub), 20)) for x in z3]
        z1n = [sp.Float(z.real, 20) + sp.I*sp.Float(z.imag, 20) for z in z1n]; z3n = [sp.Float(z.real, 20) + sp.I*sp.Float(z.imag, 20) for z in z3n]
        K = structures_direct(z1n, z3n, a)
        for s in (1, 2, 3, 4, 5):
            f = lamb(K[s])
            # y = s*yhat; kernel homogeneous of degree a-3; int_0^inf ds s^{a-1} e^{-s w} = Gamma(a) w^{-a};  d^3y = s^2 ds sinh t dt dphi
            vals = f(yh[0], yh[1], yh[2], P0a[:, None], P1a[:, None], P2a[:, None])*np.sinh(T)*wdot**(-av)
            num = _G(av)*np.sum(vals*Wt)
            Mv = M[(name, s)].subs({a: sp.Rational(9, 2), nu: sp.Rational(11, 10)}).subs(sp.Symbol('g1'), gam(av-1)/gam(av))
            Mv = complex(sp.N(Mv, 20))
            x0 = np.cosh(xiv)
            harm = Pd(x0, 0) if name.startswith('0') else np.sinh(xiv)**2*Pd(x0, 2)
            pred = 2**(av-1)*np.sqrt(np.pi)*gam(av)*Mv*harm
            print(f"   {name} D{s}:  numeric {num:.9f}   predicted {pred:.9f}")

    print("(4) stress-tensor pair, non-derivative structures: (z1,z3) = (eps_+^2, eps_-^2) [weight 0] and (eps_+^2, eps_+^2) [weight 4]")
    g1 = sp.Symbol('g1')
    def TT(z1, z3):
        # D1^2 Phi_a = (z1.z3)^2 Phi_a + 2 (z1.z3)(z1.y)(z3.y) Phi_{a-2} + (z1.y)^2 (z3.y)^2 Phi_{a-4};  D1 D2, D2^2 likewise
        r2 = sp.Rational(1, 4)*4/((a-sp.Rational(5, 2))**2 + nu**2)                 # C gamma ratio for b = a-2
        r4 = r2*sp.Rational(1, 4)*4/((a-sp.Rational(9, 2))**2 + nu**2)              # b = a-4
        T0, T2, T4 = Psi(a), r2*Psi(a-2), r4*Psi(a-4)
        q2 = zD(z1, zD(z3, T2)); q4 = zD(z1, zD(z1, zD(z3, zD(z3, T4))))
        zz = dot(z1, z3)
        return {'D1^2': zz**2*T0 + 2*zz*q2 + q4, 'D1 D2': -zz*q2 - q4, 'D2^2': q4}
    for name, (z1, z3) in {'0': (EP, EM), '4': (EP, EP)}.items():
        for s, e in TT(z1, z3).items():
            e, P0, P1 = at_point(e)
            if name == '0':
                assert sp.cancel(sp.diff(e, P1)) == 0
                print(f"   weight 0, {s}:  M = {sp.factor(sp.cancel(sp.diff(e, P0)))}")
            else:
                # P^{(4)} := sinh^4 xi d^4P/dc^4, expressed through P, P'
                P4, _, _ = reduce_P(SH**4*sp.Derivative(Pf(c), (c, 4)), CH)
                P4 = sp.cancel(P4)
                Mc = sp.cancel(sp.diff(e, P1)/sp.diff(P4, P1))
                assert sp.cancel(e - Mc*P4) == 0
                print(f"   weight 4, {s}:  M = {sp.factor(Mc)}   (times P^(4) = sinh^4 xi P'''')")
