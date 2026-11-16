"""
Caron-Huot--Li helicity structures in momentum space (arXiv:2102.08160, eq. Tpm and the E/O normalisation eq. 'even odd'), compared with
the Fourier transform of differential-basis structures (ft_helicity), to obtain N_h(Delta,J) in  T^{E,h} = N_h * e_h  exactly.
Frames: CHL put O_1(0), O_2(x), O_3(inf) and transform in x;  ours (infframe) O_1(x), O_2(0):  x_CHL = -x, so  FT_CHL(p) = FT_ours(-p).
   T^{h1,h2}(p) = (4pi)^{3/2} (-i sqrt2)^{J1+J2+J3} / 2^{tau1+tau2-Delta3} (e1 Pi_{-h1 p} e3)^{J1} (e2 Pi_{h2 p} e3)^{J2} (p.e3)^{J3-J1-J2} |p|^{beta-3},
   e Pi_{+-p} e' = (e.e' - (p.e)(p.e')/p^2 +- (i/|p|) eps(e, p, e'))/2,   beta = (Delta1+J1) + (Delta2+J2) - (Delta3+J3),
   T^{E,same} = G_E (T^{++} + T^{--})/sqrt2,  T^{E,opp} = G_E (T^{+-} + T^{-+})/sqrt2,  G_E = Gamma((3-tau1-tau2+Delta3+J3)/2)/Gamma((tau1+tau2-tau3)/2).
Our structure:  F = (-zeta.x/2)^J |x|^{-2al} R;  its FT, full normalisation:  FT = sum terms  K(e,t) ... (ft_helicity.ft_term with K including Gamma's).
"""
import sympy as sp
import mpmath as mp
from ft_helicity import R_inv, AX, BX, WX, XX, AB, AW, BW, AP, BP, WP, PP, _deriv
mp.mp.dps = 40

def eps3(a, b, c):
    return sp.Matrix.hstack(a, b, c).det()

def chl_components(l, Dv, Jv):
    """T^{E,h}(p) components at p = (0,0,-1) [our p -> -p], e1, e2 in {(1, +-i, 0)}, e3 = zeta = (1, 0, i).  Returns dict (ha,hb)->(same, opp)."""
    J1 = J2 = l; D1 = D2 = l + 1; tau = 1
    p = sp.Matrix([0, 0, 1]); e3 = sp.Matrix([1, 0, sp.I])          # CHL momentum (ours is evaluated at -p)
    d = lambda u, v: (u.T*v)[0]
    pnorm = 1
    def proj(e, s):     # e Pi_{s p} e3
        return (d(e, e3) - d(p, e)*d(p, e3)/pnorm**2 + s*sp.I/pnorm*eps3(e, p, e3))/2
    beta = (D1 + J1) + (D2 + J2) - (Dv + Jv)
    pref = (4*sp.pi)**sp.Rational(3, 2)*(-sp.I*sp.sqrt(2))**(J1 + J2 + Jv)/2**(2*tau - Dv)
    GE = sp.gamma((3 - 2*tau + Dv + Jv)/2)/sp.gamma((2*tau - (Dv - Jv))/2)
    eps = {1: sp.Matrix([1, sp.I, 0]), -1: sp.Matrix([1, -sp.I, 0])}
    out = {}
    for ha in (1, -1):
        for hb in (1, -1):
            e1, e2 = eps[ha], eps[hb]
            T = {}
            for (h1, h2) in ((1, 1), (-1, -1), (1, -1), (-1, 1)):
                T[(h1, h2)] = pref*proj(e1, -h1)**J1*proj(e2, h2)**J2*d(p, e3)**(Jv - J1 - J2)*pnorm**(beta - 3)
            same = GE*(T[(1, 1)] + T[(-1, -1)])/sp.sqrt(2)
            opp = GE*(T[(1, -1)] + T[(-1, 1)])/sp.sqrt(2)
            out[(ha, hb)] = (complex(sp.N(same, 30)), complex(sp.N(opp, 30)))
    return out

def our_components(l, coeffs, Dv, Jv):
    """full FT of F at p = (0,0,1) [evaluated at -p via the sign of p below], same e1, e2, zeta; numeric (30 digits)."""
    R = R_inv(l, coeffs, Dv, Jv)
    al = sp.nsimplify((4*l + 2 - Dv + Jv)/sp.Integer(2))
    num, den = sp.fraction(sp.together(R))
    P = sp.Poly(den, WX, XX); (dW, dX), dcoef = P.terms()[0]
    acc = 0
    for term in sp.Add.make_args(sp.expand(num)):
        pw = term.as_powers_dict()
        i, j = int(pw.get(AX, 0)), int(pw.get(BX, 0)); kW, kX = int(pw.get(WX, 0)), int(pw.get(XX, 0))
        rest = term/(AX**i*BX**j*WX**kW*XX**kX)
        e = Jv + kW - dW; t = -al + kX - dX
        K = (-sp.I)**e*sp.pi**sp.Rational(3, 2)*2**(2*t + 3 + e)*sp.gamma(t + e + sp.Rational(3, 2))/sp.gamma(-t)
        f = K*WP**e*PP**(-t - sp.Rational(3, 2) - e)
        for _ in range(i): f = sp.I*_deriv(f, 'a')
        for _ in range(j): f = sp.I*_deriv(f, 'b')
        acc += rest/dcoef*f
    acc = acc*sp.Rational(-1, 2)**Jv
    p = sp.Matrix([0, 0, -1]); zeta = sp.Matrix([1, 0, sp.I])
    d = lambda u, v: (u.T*v)[0]
    eps = {1: sp.Matrix([1, sp.I, 0]), -1: sp.Matrix([1, -sp.I, 0])}
    out = {}
    for ha in (1, -1):
        for hb in (1, -1):
            a, b = eps[ha], eps[hb]
            val = acc.subs({AP: d(a, p), BP: d(b, p), PP: 1, AB: d(a, b), AW: d(a, zeta), BW: d(b, zeta), WP: d(zeta, p)})
            out[(ha, hb)] = complex(sp.N(val, 30))
    return out
