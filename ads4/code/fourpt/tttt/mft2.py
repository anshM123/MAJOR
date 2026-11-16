"""
Exact MFT OPE data of [O O]_{N,J} (O conserved, spin l) in the differential-basis normalisation:
structures F_c(inf) = (-zeta.x/2)^J |x|^{-2 al} sum_w c_w R_w  at Delta = 2(l+1) + N  (infframe), zeta = (1, i, 0);
G_ab = <F_a(inf), F_b(inf)>_Fock (fock2),  P = (zeta.zetabar)^J G^{-1} = 2^J G^{-1}.
"""
import sympy as sp
from infframe import X, Z1, Z2, W, dot, R_words, alpha_ref, Delta, J, words
import fock2

ZETA = (1, sp.I, 0)

def structure_inf(l, coeffs, N, Jv):
    """polynomial F(x; a, b) for the coefficient dict coeffs (values may depend on Delta, J) at Delta = 2(l+1) + N, J = Jv.
    Built in the invariants (ax, bx, wx = zeta.x, xx, ab, aw, bw), then mapped to components with zeta = (1, i, 0)."""
    from ft_helicity import R_inv, AX, BX, WX, XX, AB, AW, BW
    Dv = 2*(l + 1) + N
    R = R_inv(l, coeffs, Dv, Jv)
    al = sp.Rational(4*l + 2 - Dv + Jv, 2)
    F = sp.expand(R*sp.Rational(-1, 2)**Jv*WX**Jv*XX**(-al))
    for term in sp.Add.make_args(F):
        for b, e in term.as_powers_dict().items():
            if b in (WX, XX) and (e < 0 or not sp.sympify(e).is_integer):
                raise AssertionError("non-polynomial term %s" % term)
    zs = ZETA
    csub = {AX: dot(Z1, X), BX: dot(Z2, X), WX: dot(zs, X), XX: dot(X, X), AB: dot(Z1, Z2), AW: dot(Z1, zs), BW: dot(Z2, zs)}
    return sp.expand(F.xreplace(csub))

def gram(l, structs):
    n = len(structs)
    G = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            G[i, j] = fock2.bilinear(structs[i], structs[j], X, l)
    return G
