"""
Fourier-Laplace transform of nu-harmonic tensor kernels on H^2 (the spinning generalisation of CHMRS (3.18)/(B.2)).

Kernel modes.  For yhat in H^2 and a boundary null vector l, let n_l = l/(-yhat.l) - yhat (unit tangent at yhat),
n_l^perp_mu = eps_{mu a b} yhat^a n_l^b,  z(th) = yhat + cos th n_l + sin th n_l^perp.  The (q1,q3) mode functional of
a tensor t (polynomial t(z_1, z_3) of bidegree (l,l)) is  k_{q1 q3}[t] = (2pi)^{-2} int dth1 dth3 e^{-i(q1 th1 + q3 th3)} t(z(th1), z(th3)),
and the nu-harmonic kernel is  K^nu_{q1q3}(xhat, yhat) = (1/2pi) int D l (-xhat.l)^{-1/2+i nu} (-yhat.l)^{-1/2-i nu} k_{q1q3}(yhat, l)
   = F^nu_{q1+q3}(xi) x (mode (q1,q3) in the frame at yhat pointing to xhat).

Transform.  Master formula (derived in SPINNING_E2E.md from CHMRS (3.45)):
   int_{y>0} d^3y e^{q.y} |y|^{b-3} (-yhat.l)^{-sig} = c(b, sig) |q|^{-b} (-qhat.l)^{-sig},
   c(b, sig) = 2^{b-1} sqrt(pi) Gamma((b-sig)/2) Gamma((b+sig-1)/2)      (c(a, 1/2 + i nu) = 2^{a-1} sqrt(pi) gamma_a(nu)).
For a constant tensor  t = (z_1.eps_{h1})^l (z_3.eps_{h3})^l  with eps_pm = (n_par +- i e_perp)/sqrt2 at qhat:
   int d^3y e^{q.y} |y|^{a-3} (-yhat.l)^{-1/2-i nu} k_{q1q3}(yhat,l)[t]
        = c(a, 1/2+i nu) |q|^{-a} (-qhat.l)^{-1/2-i nu} * M^{q1q3}_{h1h3}(nu, a) * e^{i mu theta},
where l/(-qhat.l) = qhat + cos(theta) n_par + sin(theta) e_perp.  M = R + g1 * S with R, S rational in (nu, a) and
g1 = gamma_{a-1}(nu)/gamma_a(nu).   M_coeff(l, q1, q3, h1, h3) returns (R, S, mu).
"""
import sympy as sp, itertools, os, pickle
from embed import cached

nu, a = sp.symbols('nu a', real=True)
s_ = sp.Rational(1, 2) + sp.I*nu
g1 = sp.Symbol('g1')

def _poch(x, n):
    n = int(n)
    if n >= 0:
        return sp.Mul(*[(x + i) for i in range(n)]) if n > 0 else sp.Integer(1)
    return 1/sp.Mul(*[(x - i) for i in range(1, -n+1)])

def cratio(e, j):
    """c(a + e - j, s + j)/c(a, s)  as  rational * {1 or g1}"""
    X1 = (a - s_)/2; X2 = (a + s_ - 1)/2
    if e % 2 == 0:
        return sp.Integer(2)**(e - j)*_poch(X1, (e - 2*j)//2)*_poch(X2, e//2)
    return sp.Integer(2)**(e - j)*_poch(X1 - sp.Rational(1, 2), (e + 1 - 2*j)//2)*_poch(X2 - sp.Rational(1, 2), (e + 1)//2)*g1

def M_coeff(l, q1, q3, h1, h3):
    E = sp.Symbol('E')                       # e^{i theta}
    th_c = (E + 1/E)/2; th_s = (E - 1/E)/(2*sp.I)
    lv = sp.Matrix([1, th_c, th_s])          # contravariant components of l (gauge -qhat.l = 1)
    eps = {+1: sp.Matrix([0, 1, sp.I])/sp.sqrt(2), -1: sp.Matrix([0, 1, -sp.I])/sp.sqrt(2)}
    eta = sp.diag(-1, 1, 1)
    def dot(u, v): return (u.T*eta*v)[0, 0]
    # symbols for the y-dependent scalars
    R, Linv = sp.symbols('R Linv')           # |y|,  1/(-y.l)
    Ys = {h: sp.Symbol('Y%s' % ('p' if h > 0 else 'm')) for h in (1, -1)}      # y.eps_h
    Ds = {h: sp.Symbol('D%s' % ('p' if h > 0 else 'm')) for h in (1, -1)}      # det[eps_h, y, l]
    def ABC(h):
        e = eps[h]
        A = Ys[h]/R
        npar = R*dot(lv, e)*Linv - Ys[h]/R
        nperp = Ds[h]*Linv
        return A, (npar - sp.I*nperp)/2, (npar + sp.I*nperp)/2       # z(th).eps = A + B e^{i th} + C e^{-i th}
    def mode(h, q):
        A, B, C = ABC(h)
        th = sp.Symbol('th')
        poly = sp.expand((A + B*th + C/th)**l)
        return sp.expand(poly).coeff(th, q)
    integrand = sp.expand(mode(h1, q1)*mode(h3, q3))
    if integrand == 0:
        return sp.Integer(0), sp.Integer(0), 0
    # directional derivative vectors (contravariant): y.eps -> eps ; det[eps, y, l] = y^al W_al, W_al = eps_{mu al be} eps^mu l^be
    def Wvec(e):
        W = sp.Matrix([0, 0, 0])
        for al in range(3):
            W[al] = sum(sp.LeviCivita(mu, al, be)*e[mu]*lv[be] for mu in range(3) for be in range(3))
        return eta*W                          # raise the index
    dirs = {Ys[1]: eps[1], Ys[-1]: eps[-1], Ds[1]: Wvec(eps[1]), Ds[-1]: Wvec(eps[-1])}
    q = sp.Matrix(sp.symbols('q0 q1 q2'))
    kap, sig = sp.symbols('kap sig')
    qq = -dot(q, q); ql = -dot(q, lv)
    def Dv(P, v):
        # derivative of  qq^kap ql^{-sig} P  divided by qq^kap ql^{-sig}
        dP = sum(v[i]*sp.diff(P, q[i]) for i in range(3))
        dqq = -2*dot(q, v); dql = -dot(lv, v)
        return dP + P*(kap*dqq/qq - sig*dql/ql)
    shift = 2*l
    poly = sp.Poly(sp.expand(integrand*R**shift), R, Linv, Ys[1], Ys[-1], Ds[1], Ds[-1])
    acc = {}                                  # power of E -> {'r': expr, 'g': expr}
    qsub = {q[0]: 1, q[1]: 0, q[2]: 0}        # |q| = 1 (the result is homogeneous of degree -a in q)
    for (eR, j, yp, ym, dp, dm), coef in poly.terms():
        eR = eR - shift
        # monomial: coef * R^eR Linv^j * (y-linear factors);  integrand |y|^{a-3} (-yhat.l)^{-s} -> |y|^{a+eR+s-3} (-y.l)^{-s-j}
        # I(b', sigma) = c(b'-sigma, sigma) |q|^{2 sigma - b'} (-q.l)^{-sigma},  b' = a + eR + s,  sigma = s + j
        P = sp.Integer(1)
        for sym, cnt in ((Ys[1], yp), (Ys[-1], ym), (Ds[1], dp), (Ds[-1], dm)):
            for _ in range(cnt):
                P = Dv(P, dirs[sym])
        kval = (2*(s_ + j) - (a + eR + s_))/2
        P = sp.expand((coef*P).subs({kap: kval, sig: s_ + j}).subs(qsub))
        cr = cratio(eR, j)
        kind = 'g' if cr.has(g1) else 'r'
        cr = cr.subs(g1, 1)
        nshift = 4*l + 4
        pE = sp.Poly(sp.expand(P*E**nshift), E)
        for (pw,), cE in pE.terms():
            d = acc.setdefault(pw - nshift, {'r': sp.Integer(0), 'g': sp.Integer(0)})
            d[kind] += cr*cE
    res = {}
    for pw, d in acc.items():
        Rr = sp.factor(sp.cancel(sp.together(d['r']))); S = sp.factor(sp.cancel(sp.together(d['g'])))
        if Rr != 0 or S != 0:
            res[pw] = (Rr, S)
    if len(res) == 0:
        return sp.Integer(0), sp.Integer(0), 0
    assert len(res) == 1, res
    mu = list(res)[0]
    return res[mu][0], res[mu][1], mu

def M_table(l):
    def build():
        out = {}
        for q1 in range(-l, l+1):
            for q3 in range(-l, l+1):
                for h1 in (1, -1):
                    for h3 in (1, -1):
                        out[(q1, q3, h1, h3)] = M_coeff(l, q1, q3, h1, h3)
        return out
    return cached('Mtable_l%d' % l, build)

if __name__ == "__main__":
    import time
    t0 = time.time()
    for l in (1, 2):
        T = M_table(l)
        print("l =", l, " (%.0f s)" % (time.time() - t0))
        for k, v in sorted(T.items()):
            if v[0] != 0 or v[1] != 0:
                print("  (q1,q3,h1,h3) =", k, " mu =", v[2], "  R =", v[0], "  S =", v[1])
