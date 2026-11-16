"""
Symbolic (in Delta, J) proof attempt of  kappa_h A_h = 2^{3J+3/2}  (t8_kappa_exact done with symbols instead of rational points).
X_h(Delta,J) (chl_xspace, symbolic exponents) is decomposed into the differential-basis words over the field Q(Delta,J) (l = 2: plus the
Gram-determinant term mu G4 XX^2/WX^2 of the d = 3 relation); T is applied; lambda_h is extracted symbolically and the identity
   lambda_h C_l(Delta,J)/C_l(3-Delta,-1-J) * A_h = 2^{3J+3/2}
is checked as an identity of rational functions (powers of 2 separated by hand: C_l carries 2^{J/2}; the words are normalised by 2^{-J}).
Usage: python3 t10_kappa_symbolic.py l     Output: out_t10_kappa_symbolic_l<l>.txt
"""
import os, sys, time
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from chl_xspace import xspace_structure
from ft_helicity import R_inv, AX, BX, WX, XX, AB, AW, BW
from infframe import words
from embed2 import shifts, Delta as sD, J as sJ
from t1_structures_l1 import b_ratio
INV = (AX, BX, WX, XX, AB, AW, BW)
D, Jj = sp.symbols('D J')
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def strip(expr, Jv, al):
    """divide each term by WX^J XX^{-al} (symbolic exponents)"""
    res = 0
    for t in sp.Add.make_args(sp.expand(expr)):
        c, facs = t.as_coeff_mul()
        ew = 0; ex = 0; rest = c
        for f in facs:
            b, e = f.as_base_exp()
            if b == WX: ew += e
            elif b == XX: ex += e
            else: rest *= f
        res += rest*WX**sp.expand(ew - Jv)*XX**sp.expand(ex + al)
    return sp.expand(res)

def decompose(l, Dv, Jv, kind):
    al = (4*l + 2 - Dv + Jv)/2
    Xn = strip(xspace_structure(l, Dv, Jv, kind), Jv, al)
    for t in sp.Add.make_args(Xn):
        for b, e in t.as_powers_dict().items():
            if b in (WX, XX): assert e.is_integer, t
    ws = words(l)
    Rw = [sp.together(R_inv(l, {w: 1}, Dv, Jv).xreplace({AX: -AX, BX: -BX, WX: -WX})) for w in ws]   # words x 2^J
    if l == 2:
        G4 = sp.Matrix([[XX, AX, BX, WX], [AX, 0, AB, AW], [BX, AB, 0, BW], [WX, AW, BW, 0]]).det()
        Rw.append(G4*XX**2/WX**2)
    n = len(Rw)
    # common denominator: monomial in WX, XX
    allx = [Xn] + Rw
    den = 1
    for e in allx:
        d = sp.fraction(sp.together(e))[1]
        den = sp.lcm(den, d)
    cols = [sp.expand(sp.cancel(e*den)) for e in allx]
    mons = set()
    polys = []
    for c in cols:
        P = sp.Poly(c, *INV); polys.append(P.as_dict()); mons |= set(P.as_dict())
    mons = sorted(mons)
    K = sp.QQ.frac_field(D, Jj)
    A = [[K.from_sympy(sp.sympify(polys[j + 1].get(m, 0))) for j in range(n)] for m in mons]
    b = [[K.from_sympy(sp.sympify(polys[0].get(m, 0)))] for m in mons]
    Am = DomainMatrix(A, (len(mons), n), K); bm = DomainMatrix(b, (len(mons), 1), K)
    aug = Am.hstack(bm)
    R, piv = aug.rref()
    assert n not in piv, "inconsistent"
    assert len(piv) == n, "word coefficients not unique (%d pivots)" % len(piv)
    sol = [K.to_sympy(R[i, n].element) for i in range(n)]
    return ws, sol[:len(ws)]

if __name__ == "__main__":
    l = int(sys.argv[1])
    for kind in ('same', 'opp'):
        t = time.time()
        ws, c1 = decompose(l, D, Jj, kind)
        pr(f"  l={l} {kind}: decomposition at (Delta,J) done ({time.time()-t:.0f} s)")
        Dt, Jt = 3 - D, -1 - Jj
        _, c2 = decompose(l, Dt, Jt, kind)
        pr(f"  l={l} {kind}: decomposition at (3-Delta,-1-J) done ({time.time()-t:.0f} s)")
        Tc = [sp.cancel(c1[i]/b_ratio(shifts(w)[0] - shifts(w)[1]).subs({sD: D, sJ: Jj}, simultaneous=True)) for i, w in enumerate(ws)]
        i0 = max(range(len(ws)), key=lambda i: len(str(c2[i])) if c2[i] != 0 else -1)
        lam = sp.cancel(Tc[i0]/c2[i0])                       # lambda_h / 2^{J - J~} (both decompositions carry 2^J, 2^J~)
        resid = [sp.cancel(Tc[i] - lam*c2[i]) for i in range(len(ws))]
        ok = all(r == 0 for r in resid)
        # kappa A / 2^{3J+3/2} = lam * 2^{2J+1} * C(D)/C(Dt) * A / 2^{3J+3/2};  C ratio = 2^{J+1/2} * Pochhammer ratio
        rf = lambda x, n: sp.Mul(*[x + i for i in range(n)])
        Cp = lambda Dv, Jv: rf((2 - Dv + Jv)/2, 2*l)*rf((Dv + Jv + 1 - 4*l)/2, 2*l)
        A = rf(-Jj, 2*l)/rf(Jj + 1, 2*l) if kind == 'opp' else rf(D - 1, 2*l)/rf(2 - D, 2*l)
        test = sp.factor(sp.cancel(lam*Cp(D, Jj)/Cp(Dt, Jt)*A))      # powers of 2: 2^{2J+1} 2^{J+1/2} = 2^{3J+3/2}
        pr(f"  l={l} {kind}: T X_h(Delta) = lambda X_h(Delta~) identically: {ok};  lambda/2^(2J+1) = {sp.factor(lam)};  "
           f"kappa_h A_h / 2^(3J+3/2) = {test}   ({time.time()-t:.0f} s)")
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_t10_kappa_symbolic_l%d.txt' % l), 'w').write("\n".join(out) + "\n")
