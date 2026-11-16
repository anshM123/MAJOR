"""
Caron-Huot--Li parity-even helicity structures T^{E,h} (arXiv:2102.08160, eqs. Tpm, even odd) in POSITION space, for two conserved spin-l
operators, by term-wise inverse Fourier transform of the momentum-space definition; exact (Gamma functions reduced to Pochhammer ratios,
one common scalar factor C(Delta,J) dropped).

Momentum space (CHL frame O_1(0), O_2(x), O_3(inf), FT in x):  with A_i = e_i.e3 - (p.e_i)(p.e3)/p^2, B_i = eps(e_i, p, e3),
   e_i Pi_{s p} e3 = (A_i + s i B_i/|p|)/2,
   T^{h1,h2} = pref (e1 Pi_{-h1 p} e3)^l (e2 Pi_{h2 p} e3)^l (p.e3)^{J-2l} |p|^{beta-3},   beta = 2(2l+1) - (Delta + J)  [Delta_i = l+1, J_i = l],
   same = T^{++} + T^{--},  opp = T^{+-} + T^{-+}  (common factors G_E/sqrt2, pref dropped: they are h-independent).
Products of B's: B1 B2, B1^2, B2^2 are Gram determinants (e_i null).
Inverse FT:  int d^3p/(2pi)^3 e^{ipx} (e3.p)^n |p|^{2s} = i^n (2pi)^{-3} pi^{3/2} 2^{2s+3+n} Gamma(s+n+3/2)/Gamma(-s) (e3.x)^n |x|^{-2s-3-2n},
   (p.e_i) -> -i e_i.d_x.
Output: the structure in the invariants  AX = e1.x, BX = e2.x, WX = e3.x, XX = x^2, AB = e1.e2, AW = e1.e3, BW = e2.e3  of the CHL frame;
our frame (infframe: O_1(x), O_2(0)) is x_ours = -x_CHL.
"""
import sympy as sp
from ft_helicity import AX, BX, WX, XX, AB, AW, BW

P1, P2, P3, PP = sp.symbols('P1 P2 P3 PP')        # p.e1, p.e2, p.e3, p^2

def _gr(x, n):
    n = int(n)
    return sp.rf(x, n) if n >= 0 else 1/sp.rf(x + n, -n)

def momentum_structure(l, Dv, Jv, kind):
    """polynomial in P1, P2, P3 (P3 also to the power J - 2l, kept symbolic as an overall P3^J factor) and PP^s; returns dict
    (a, b, k, s_shift) -> coefficient for terms  P1^a P2^b P3^(J + k) PP^(s0 + s_shift),  s0 = (beta - 3)/2."""
    # A_i, B-products in the invariants
    A1 = AW - P1*P3/PP; A2 = BW - P2*P3/PP
    G12 = sp.Matrix([[AB, P1, AW], [P2, PP, P3], [BW, P3, 0]]).det()          # B1 B2
    G11 = sp.Matrix([[0, P1, AW], [P1, PP, P3], [AW, P3, 0]]).det()           # B1^2
    G22 = sp.Matrix([[0, P2, BW], [P2, PP, P3], [BW, P3, 0]]).det()           # B2^2
    b1, b2 = sp.symbols('b1 b2')
    combos = [(-1, 1), (1, -1)] if kind == 'same' else [(-1, -1), (1, 1)]      # (sign for e1, sign for e2) in Pi
    tot = 0
    for (s1, s2) in combos:
        tot += sp.expand((A1 + s1*sp.I*b1)**l*(A2 + s2*sp.I*b2)**l)
    # b_i = B_i/|p| ; reduce products
    poly = sp.Poly(tot, b1, b2)
    red = 0
    for (i, j), c in poly.terms():
        assert (i + j) % 2 == 0
        if i % 2 == 0 and j % 2 == 0:
            red += c*G11**(i//2)*G22**(j//2)/PP**((i + j)//2)
        else:
            red += c*G11**((i-1)//2)*G22**((j-1)//2)*G12/PP**((i + j)//2)
    red = sp.expand(red*P3**(-2*l))           # times (p.e3)^{J-2l}: P3^{-2l} here, P3^J overall
    out = {}
    for term in sp.Add.make_args(red):
        pw = term.as_powers_dict()
        a, b, k, s = int(pw.get(P1, 0)), int(pw.get(P2, 0)), int(pw.get(P3, 0)), int(pw.get(PP, 0))
        c = term/(P1**a*P2**b*P3**k*PP**s)
        out[(a, b, k, s)] = out.get((a, b, k, s), 0) + c
    return out

def _deriv_x(expr, which):
    """e1.d_x (which = 'a') or e2.d_x ('b') on a function of AX, BX, WX, XX (e1, e2 null)"""
    if which == 'a':
        rules = {AX: 0, BX: AB, WX: AW, XX: 2*AX}
    else:
        rules = {AX: AB, BX: 0, WX: BW, XX: 2*BX}
    return sum(sp.diff(expr, v)*r for v, r in rules.items() if r != 0)

def xspace_structure(l, Dv, Jv, kind):
    """x-space T^{E,kind} (CHL frame) up to a common h-independent scalar, as an expression in AX, BX, WX, XX, AB, AW, BW
    (exact; Delta, J rational).  Common factor dropped: i^J (2pi)^{-3} pi^{3/2} 2^{2 s0 + 3 + J} Gamma(s0 + J + 3/2)/Gamma(-s0) * G_E pref/sqrt2."""
    beta = 2*(2*l + 1) - (Dv + Jv)
    s0 = sp.nsimplify((beta - 3)/sp.Integer(2))
    terms = momentum_structure(l, Dv, Jv, kind)
    acc = 0
    for (a, b, k, s), c in terms.items():
        n = Jv + k; ss = s0 + s
        dn, ds = k, s
        # ratio of the inverse-FT constant to the reference (n0 = J, s0): i^{dn} 2^{2ds + dn} Gamma(ss+n+3/2)/Gamma(s0+J+3/2) * Gamma(-s0)/Gamma(-ss)
        Kr = sp.I**dn*sp.Integer(2)**(2*ds + dn)*_gr(s0 + Jv + sp.Rational(3, 2), ds + dn)/_gr(-s0, -ds)
        f = Kr*WX**n*XX**(-ss - sp.Rational(3, 2) - n)
        for _ in range(a): f = -sp.I*_deriv_x(f, 'a')
        for _ in range(b): f = -sp.I*_deriv_x(f, 'b')
        acc += c*f
    return sp.expand(acc)

def to_our_frame(expr):
    """x_ours = -x_CHL: AX, BX, WX -> -AX, -BX, -WX"""
    return expr.xreplace({AX: -AX, BX: -BX, WX: -WX})

EPS = sp.Symbol('eps')

def _lead(r):
    """leading behaviour c eps^k of a rational function of EPS at 0: returns (k, c)"""
    n, d = sp.fraction(sp.cancel(sp.together(r)))
    pn = sp.Poly(n, EPS); pd = sp.Poly(d, EPS)
    kn = min(m[0] for m in pn.monoms()); kd = min(m[0] for m in pd.monoms())
    return kn - kd, pn.coeff_monomial(EPS**kn)/pd.coeff_monomial(EPS**kd)

def xspace_structure_limit(l, Dv, Jv, kind):
    """as xspace_structure, but at Delta -> Dv + eps' (eps = -eps'/2 shift of s0) and returning the leading term in eps:
    (k, S) with T^{E,kind}(Dv + eps') = C(Delta) eps^k (S + O(eps)) [common XX^{-eps} dropped].  Needed at double-twist points where
    individual Gamma ratios have poles (k < 0) or zeros (contact terms)."""
    beta = 2*(2*l + 1) - (Dv + Jv)
    s0n = sp.nsimplify((beta - 3)/sp.Integer(2))
    s0 = s0n + EPS
    terms = momentum_structure(l, Dv, Jv, kind)
    lead = []
    for (a, b, k, s), c in terms.items():
        n = Jv + k; ss = s0n + s
        Kr = sp.I**k*sp.Integer(2)**(2*s + k)*_gr(s0 + Jv + sp.Rational(3, 2), s + k)/_gr(-s0, -s)
        kk, cc = _lead(Kr)
        f = cc*WX**n*XX**(-ss - sp.Rational(3, 2) - n)
        for _ in range(a): f = -sp.I*_deriv_x(f, 'a')
        for _ in range(b): f = -sp.I*_deriv_x(f, 'b')
        lead.append((kk, c*f))
    korders = sorted(set(k for k, _ in lead))
    for k0 in korders:
        S = sp.expand(sum(f for k, f in lead if k == k0))
        if S != 0:
            return k0, S
    return None, 0
