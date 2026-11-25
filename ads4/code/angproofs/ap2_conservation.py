"""Step 2 (proof): the candidate symbolic-J structures of ap1_fit are conserved at x1 and Bose symmetric under
(x1,z1) <-> (x3,z3), as exact rational-function identities in (x1, z1, a3, a4, s, J).

Frame: x4 = 0, x3 = (0,0,1) (any pair of spacelike-separated points is conformally equivalent to this one, and the
conservation operator at x1 is conformally covariant because Delta_T = 3 = J_T + d - 2), z3 = zeta(a3,1), z4 = zeta(a4,1)
(generic null polarisations up to scale; the structures are polynomial in z3, z4), x1 and z1 fully symbolic (z1 off the
cone: the Todorov operator is interior).  The power prefactor V4^J X13^{-e13} X14^{-e14} X34^{-e34}, with SYMBOLIC
exponents, is differentiated analytically (logarithmic derivative); all remaining exponents are integers.
usage: python3 ap2_conservation.py [name ...]   (pp pm odd oddJ)"""
import pickle, sys, time
from common import *
import structures as st

x1 = list(sp.symbols('x1_0:3')); z1 = list(sp.symbols('z1_0:3'))
a3, a4 = sp.symbols('a3 a4')
ETA = (-1, 1, 1)
def zeta(a, b):
    return [a*b, (a*a + b*b)/2, sp.I*(a*a - b*b)/2]
X3 = [sp.Integer(0), sp.Integer(0), sp.Integer(1)]; X4 = [sp.Integer(0)]*3
Z3, Z4 = zeta(a3, sp.Integer(1)), zeta(a4, sp.Integer(1))

def R_of(mr, B, Jv=None):
    m1, m3, k, n13, n14, n34, eps = mr
    r = B['V1']**m1*B['V3']**m3*B['V4']**(-k)*B['H13']**n13*B['H14']**n14*B['H34']**n34
    if eps:
        r = r*B['eps' + eps]*B['X' + eps]
    return r

def check_conservation(vec, parity):
    half = sp.Rational(1, 2) if parity == 'odd' else 0
    e13, e14 = 5 - S + half, S + half
    B = st.blocks(x1, z1, X3, Z3, X4, Z4)
    B = {k: sp.together(v) for k, v in B.items()}
    dlog = [J*sp.diff(B['V4'], x1[m])/B['V4'] - e13*sp.diff(B['X13'], x1[m])/B['X13']
            - e14*sp.diff(B['X14'], x1[m])/B['X14'] for m in range(3)]
    R = sum(c*R_of(mr, B) for mr, c in vec.items() if c != 0)
    dR = [sp.diff(R, z1[m]) for m in range(3)]
    box = sum(ETA[m]*sp.diff(dR[m], z1[m]) for m in range(3))
    DR = [sp.Rational(3, 2)*ETA[m]*dR[m] - sp.Rational(1, 2)*z1[m]*box for m in range(3)]
    div = sum(sp.diff(DR[m], x1[m]) + dlog[m]*DR[m] for m in range(3))
    # div is linear in z1; check each coefficient vanishes identically
    out = []
    for nu in range(3):
        c = sp.diff(div, z1[nu]).subs({z: 0 for z in z1})
        num = sp.numer(sp.together(c))
        out.append(sp.expand(num) == 0)
    return out

def check_symmetry(vec, parity, Jpar):
    """S(1,3) = S(3,1); with x1 symbolic, x3 = (0,0,1), x4 = 0, and z1 = zeta(a1,1)."""
    a1 = sp.Symbol('a1'); Z1 = zeta(a1, sp.Integer(1))
    B = st.blocks(x1, Z1, X3, Z3, X4, Z4)
    Bs = st.blocks(X3, Z3, x1, Z1, X4, Z4)
    assert sp.simplify(Bs['V4'] + B['V4']) == 0 and sp.simplify(Bs['X14'] - B['X34']) == 0
    # prefactor V4^J X13^-e13 (X14 X34)^-e14 is symmetric up to (-1)^J from V4 -> -V4
    sgn = (-1)**Jpar
    tot = sum(c*(R_of(mr, B) - sgn*R_of(mr, Bs)) for mr, c in vec.items() if c != 0)
    return sp.expand(sp.numer(sp.together(tot))) == 0

if __name__ == '__main__':
    structs = pickle.load(open(os.path.join(CACHE, 'symbolic_structures.pkl'), 'rb'))
    names = sys.argv[1:] or ['pp', 'pm', 'odd', 'oddJ']
    for name in names:
        par = 'even' if name in ('pp', 'pm') else 'odd'
        t0 = time.time()
        cons = check_conservation(structs[name], par)
        t1 = time.time()
        sym = check_symmetry(structs[name], par, 1 if name == 'oddJ' else 0)
        print(f"{name:5s}: conservation (3 components) identically zero: {cons} [{t1-t0:.0f}s];  "
              f"Bose symmetry: {sym} [{time.time()-t1:.0f}s]", flush=True)
