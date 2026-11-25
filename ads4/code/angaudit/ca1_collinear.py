"""
Audit of ANGULAR_ROUTE §21 item 1: collinear all-spin identities.

Definition used (repo conventions, KERNELS_D3 §3.1): evaluate the canonical conserved structures
<O(x4, zeta4) T(x1, z1) T(x3, z3)> at three COLLINEAR points in the Euclidean spatial plane t = 0,
   x1 = (0,-1,0), x3 = (0,1,0), x4 = (0,0,0)
with O polarised by e_+ = (0,1,i) or e_- = (0,1,-i) (its top planar-spin component, a 2d quasi-primary of weights
(h,s) or (s,h) in (w, wbar), w = x^1 + i x^2), and T polarisations z = zeta(a,b) = a^2 e_+/2 + b^2 e_-/2 + ab e_0, whose
component a^{2+m} b^{2-m} is the planar-spin-m component.  The 2d three-point structure at this configuration is
2^{Delta-6} (|u13| = 2, |u14| = |u34| = 1), so F_{m1,m3} := 2^{6-Delta} * <...>|_{component}.

Prints, for J = 4..12, ratios of F to the claimed right-hand sides (symbolic in s).
"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'kernels'))
os.environ.setdefault('KERNELS_CACHE', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache'))
import sympy as sp
import structures as st

S = sp.Symbol('s')
a1, b1, a3, b3 = sp.symbols('a1 b1 a3 b3')
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'kernels', 'data')

def zeta(a, b):
    return [a*b, (a*a + b*b)/2, sp.I*(a*a - b*b)/2]

def load(J):
    d = json.load(open(os.path.join(DATA, f'structures_J{J}.json')))
    out = {}
    for name, v in d.items():
        monos = [tuple(m[:6]) + (m[6],) for m in v['monomials']]
        coeffs = [sp.sympify(c, locals={'s': S}) for c in v['coefficients']]
        out[name] = (v['parity'], monos, coeffs)
    return out

def collinear_F(J, struct, opol='+', pts=None):
    """returns dict (m1,m3) -> F (sympy expression in s)"""
    parity, monos, coeffs = struct
    x1, x3, x4 = pts or ([0, -1, 0], [0, 1, 0], [0, 0, 0])
    z4 = [0, 1, sp.I] if opol == '+' else [0, 1, -sp.I]
    B = st.blocks([sp.Integer(c) for c in x1], zeta(a1, b1), [sp.Integer(c) for c in x3], zeta(a3, b3),
                  [sp.Integer(c) for c in x4], z4)
    tot = 0
    for m, c in zip(monos, coeffs):
        if c == 0:
            continue
        tot += c*sp.expand(st.mono_value(m, B))
    tot = sp.expand(tot)
    X13, X14, X34 = B['X13'], B['X14'], B['X34']
    half = sp.Rational(1, 2) if parity == 'odd' else 0
    Delta = 2*S - J
    pref = X13**(-(5 - S + half))*X14**(-(S + half))*X34**(-(S + half))
    pref = sp.powsimp(pref*2**(6 - Delta), force=True)
    P = sp.Poly(tot, a1, b1, a3, b3)
    out = {}
    for (e1, f1, e3, f3), c in P.terms():
        m1, m3 = e1 - 2, e3 - 2
        out[(m1, m3)] = sp.simplify(sp.factor(c)*pref)
    return out

if __name__ == '__main__':
    Js = [int(x) for x in sys.argv[1:]] or [4, 5, 6, 7, 8, 9, 10, 11, 12]
    lpp, lpm, lod = sp.symbols('lpp lpm lodd')
    for J in Js:
        d = load(J)
        h = S - J
        p = (h - 1)*(h - 2)*(h - 3)
        print(f"===== J = {J}: structures {list(d)}")
        for opol in ('+', '-'):
            Ftot = {}
            for name, struct in d.items():
                lam = {'pp': lpp, 'pm': lpm, 'odd': lod, 'oddJ': lod}[name]
                for k, v in collinear_F(J, struct, opol).items():
                    Ftot[k] = Ftot.get(k, 0) + lam*v
            for k in sorted(Ftot):
                Ftot[k] = sp.factor(sp.simplify(Ftot[k]))
            print(f"  O polarised e_{opol}:")
            for k in [(-2, -2), (-1, -1), (-2, -1), (-1, -2), (2, 2), (1, 1), (2, 1)]:
                v = Ftot.get(k, 0)
                print(f"    F{k} = {v}")
            if J % 2 == 0:
                c22 = sp.factor(sp.simplify(Ftot.get((-2, -2), 0)/(p*(h - 4))))
                c11 = sp.factor(sp.simplify(Ftot.get((-1, -1), 0)/(-16*p)))
                c21 = sp.factor(sp.simplify(Ftot.get((-2, -1), 0)/(-2*(2*h - 7)*(2*h - 5)*(2*h - 3))))
                print(f"    F(-2,-2)/[p(h)(h-4)]     = {c22}")
                print(f"    F(-1,-1)/[-16 p(h)]      = {c11}")
                print(f"    F(-2,-1)/[-2(2h-7)(2h-5)(2h-3)] = {c21}")
