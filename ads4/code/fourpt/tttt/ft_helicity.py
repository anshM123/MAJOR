"""
Helicity content of conserved <O_1(x) O_2(0) X(inf, zeta)> structures in momentum space (the definition of Caron-Huot--Li,
arXiv:2102.08160, eq. Tpm: Fourier transform in the position of one operator, helicity projections transverse to p).

Structure (infframe):  F = (-zeta.x/2)^J |x|^{-2 al} R(x; a, b; zeta),  R polynomial/rational in the invariants
   ax = a.x, bx = b.x, wx = zeta.x, xx = x^2, ab, aw, bw   (a, b, zeta null).
Term-wise Fourier transform (d = 3, Euclidean, analytic continuation in the exponents):
   FT[(zeta.x)^e |x|^{2t}] = (-i)^e pi^{3/2} 2^{2t+3+e} Gamma(t+e+3/2)/Gamma(-t) (zeta.p)^e |p|^{-2t-3-2e}
   (harmonic polynomial (zeta.x)^e, zeta null; Bochner), and (a.x) -> i a.d_p, (b.x) -> i b.d_p.
Evaluated at p = (0,0,1) with a, b in {eps_+, eps_-} (transverse, null) -> the 4 helicity components T(h_a, h_b).
"""
import sympy as sp
from infframe import word_exprs, Delta, J, words

AX, BX, WX, XX, AB, AW, BW = sp.symbols('AX BX WX XX AB AW BW')
AP, BP, WP, PP = sp.symbols('AP BP WP PP')
ISUB = None

def inv_sub(S):
    return {S['s12']: -XX, S['s15']: -1, S['s25']: -1, S['z12']: -AX, S['z21']: BX, S['z15']: 0, S['z25']: 0,
            S['z51']: WX, S['z52']: 0, S['z1z2']: AB, S['z1z5']: AW, S['z2z5']: BW}

def R_inv(l, coeffs, Dv, Jv):
    ex, S = word_exprs(l)
    sub = inv_sub(S)
    num = {Delta: Dv, J: Jv}
    tot = sum(sp.sympify(c).subs(num)*ex[w].subs(num).subs(sub) for w, c in coeffs.items())
    return sp.expand(sp.together(tot))

def _deriv(expr, which):
    """a.d_p or b.d_p acting on a function of AP, BP, WP, PP (a, b null)"""
    if which == 'a':
        rules = {AP: 0, BP: AB, WP: AW, PP: 2*AP}
    else:
        rules = {AP: AB, BP: 0, WP: BW, PP: 2*BP}
    return sum(sp.diff(expr, s)*r for s, r in rules.items() if r != 0)

def ft_term(i, j, e, t):
    """FT of (a.x)^i (b.x)^j (zeta.x)^e (x^2)^t, up to the common factor pi^{3/2}: returns expression in AP, BP, WP, PP"""
    K = (-sp.I)**e*2**(2*t + 3 + e)*sp.gamma(t + e + sp.Rational(3, 2))/sp.gamma(-t)
    f = K*WP**e*PP**(-t - sp.Rational(3, 2) - e)
    for _ in range(i):
        f = sp.I*_deriv(f, 'a')
    for _ in range(j):
        f = sp.I*_deriv(f, 'b')
    return f

def helicity_components(l, coeffs, Dv, Jv):
    """T(h_a, h_b) for h in {+1, -1} at p = (0,0,1), zeta = (1, 0, i) [zeta.p = i]."""
    R = R_inv(l, coeffs, Dv, Jv)
    al = (4*l + 2 - Dv + Jv)/sp.Integer(2)
    # F = (-1/2)^J WX^J XX^{-al} R ;  expand R into monomials AX^i BX^j WX^k XX^m (other invariants are constants)
    R = sp.expand(R*WX**0)
    num, den = sp.fraction(sp.together(R))
    # den is a monomial in WX, XX
    terms = sp.Add.make_args(sp.expand(num))
    out = {}
    eps = {1: sp.Matrix([1, sp.I, 0])/sp.sqrt(2), -1: sp.Matrix([1, -sp.I, 0])/sp.sqrt(2)}
    zeta = sp.Matrix([1, 0, sp.I]); p = sp.Matrix([0, 0, 1])
    d = lambda u, v: (u.T*v)[0]
    # collect exponents
    dpw = sp.Poly(den, WX, XX).monoms()
    assert len(sp.Poly(den, WX, XX).terms()) == 1
    dW, dX = dpw[0]; dcoef = sp.Poly(den, WX, XX).coeffs()[0]
    acc = 0
    for term in terms:
        c, mon = term.as_coeff_mul()
        powers = {AX: 0, BX: 0, WX: 0, XX: 0}
        rest = c
        for f in mon:
            b, ex = f.as_base_exp()
            if b in powers: powers[b] += ex
            else: rest *= f
        e = Jv + powers[WX] - dW
        t = -al + powers[XX] - dX
        acc += rest/dcoef*ft_term(powers[AX], powers[BX], e, t)
    acc = acc*sp.Rational(-1, 2)**Jv
    for ha in (1, -1):
        for hb in (1, -1):
            a = eps[ha]; b = eps[hb]
            val = acc.subs({AP: d(a, p), BP: d(b, p), WP: d(zeta, p), PP: 1, AB: d(a, b), AW: d(a, zeta), BW: d(b, zeta)})
            out[(ha, hb)] = complex(sp.N(val, 30))
    return out
