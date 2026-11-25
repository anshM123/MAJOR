"""Step 3 (proof, given steps 1-2): the collinear identities for SYMBOLIC J.
Configuration (ANGULAR_AUDIT ca1): x1 = (0,-1,0), x3 = (0,1,0), x4 = 0 on a line in the Euclidean plane t = 0;
O polarised by e_+ = (0,1,i); z_i = zeta(a_i,b_i).  F_{m1,m3} := 2^{6-Delta} x [a1^{2+m1} b1^{2-m1} a3^{2+m3} b3^{2-m3}] <O T T>.
Every monomial is evaluated with V4^{J-k} = V4^J V4^{-k}; V4 = -1/2 at this configuration, so all J-dependence is
explicit (polynomial coefficients of ap1/ap2 and the factor (2 V4)^J = (-1)^J)."""
import pickle
from common import *
import structures as st

a1, b1, a3, b3 = sp.symbols('a1 b1 a3 b3')
def zeta(a, b):
    return [a*b, (a*a + b*b)/2, sp.I*(a*a - b*b)/2]
P1, P3, P4 = [0, -1, 0], [0, 1, 0], [0, 0, 0]
B = st.blocks([sp.Integer(c) for c in P1], zeta(a1, b1), [sp.Integer(c) for c in P3], zeta(a3, b3),
              [sp.Integer(c) for c in P4], [0, 1, sp.I])
B = {k: sp.factor(v) for k, v in B.items()}

def collinear(vec, parity):
    v4 = B['V4']
    half = sp.Rational(1, 2) if parity == 'odd' else 0
    tot = 0
    for (m1, m3, k, n13, n14, n34, eps), c in vec.items():
        if c == 0:
            continue
        r = B['V1']**m1*B['V3']**m3*v4**(-k)*B['H13']**n13*B['H14']**n14*B['H34']**n34
        if eps:
            r = r*B['eps' + eps]*B['X' + eps]
        tot += c*sp.expand(r)
    tot = sp.expand(tot)
    # prefactor: 2^{6-Delta} v4^J X13^{-(5-s+half)} X14^{-(s+half)} X34^{-(s+half)}, X13 = 4, X14 = X34 = 1,
    # Delta = 2s - J  =>  2^{J - 4 - 2 half} v4^J = 2^{-4-2half} (2 v4)^J
    assert B['X13'] == 4 and B['X14'] == 1 and B['X34'] == 1
    pref = sp.Integer(2)**(-4 - 2*half)
    P = sp.Poly(tot, a1, b1, a3, b3)
    return {(e1 - 2, e3 - 2): sp.factor(c*pref) for (e1, f1, e3, f3), c in P.terms()}, 2*v4

if __name__ == '__main__':
    print("building blocks at the collinear point:")
    for k in ('V1', 'V3', 'V4', 'H13', 'H14', 'H34', 'eps13', 'eps14', 'eps34', 'X13'):
        print(f"   {k:5s} = {B[k]}")
    structs = pickle.load(open(os.path.join(CACHE, 'symbolic_structures.pkl'), 'rb'))
    h = S - J
    p = (h - 1)*(h - 2)*(h - 3)
    lpp, lpm, lod = sp.symbols('lambda_pp lambda_pm lambda_odd')
    F = {}
    for name, lam in (('pp', lpp), ('pm', lpm), ('odd', lod)):
        Fn, sgn = collinear(structs[name], 'even' if name != 'odd' else 'odd')
        print(f"{name}: (2V4)^J = ({sgn})^J; nonzero components {sorted(Fn)}")
        for kk, v in Fn.items():
            F[kk] = F.get(kk, 0) + lam*v
    # even J: (2V4)^J = 1
    claims = {(-2, -2): p*(h - 4)*(lpm - 2*lpp),
              (-1, -1): -16*p*(2*(S - 1)*lpp + S*lpm),
              (-2, -1): -2*(2*h - 7)*(2*h - 5)*(2*h - 3)*lod}
    for kk, cl in claims.items():
        r = sp.simplify(F[kk]/cl)
        print(f"even J: F{kk} = {sp.factor(F[kk])}\n        ratio to the claim = {r}")
    Fo, sgn = collinear(structs['oddJ'], 'odd')
    r = sp.simplify(Fo[(-2, -1)]/(-2*(2*h - 7)*(2*h - 5)*(2*h - 3)))
    print(f"odd J: (2V4)^J = ({sgn})^J = -1;  F(-2,-1) (without the (-1)^J) = {sp.factor(Fo[(-2,-1)])};  ratio = {r} (times (-1)^J = -1)")
    print("odd J, no oddJ contribution to (-2,-2),(-1,-1):", Fo.get((-2, -2), 0), Fo.get((-1, -1), 0))
    print("odd structure (even J) in (-2,-2),(-1,-1):", collinear(structs['odd'], 'odd')[0].get((-2, -2), 0),
          collinear(structs['odd'], 'odd')[0].get((-1, -1), 0))
    print("even structures in (-2,-1):", [collinear(structs[n], 'even')[0].get((-2, -1), 0) for n in ('pp', 'pm')])
    # contributing monomials
    for name in ('pp', 'pm', 'odd', 'oddJ'):
        contrib = []
        for mr, c in structs[name].items():
            if c == 0: continue
            Fm, _ = collinear({mr: c}, 'even' if name in ('pp', 'pm') else 'odd')
            if any(Fm.get(kk, 0) != 0 for kk in claims):
                contrib.append((mr, sp.factor(c), {kk: Fm.get(kk, 0) for kk in claims if Fm.get(kk, 0) != 0}))
        print(f"monomials of {name} contributing to the three components:")
        for x in contrib:
            print("    ", x)
