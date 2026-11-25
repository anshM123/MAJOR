"""Step 2 (proof): the candidate symbolic-J structures of ap1_fit are conserved at x1 and Bose symmetric under
(x1,z1) <-> (x3,z3), as exact identities in (z1, a3, a4, s, J).

Reduction to one configuration.  The divergence  Div = d_{x1}.D_{z1} S  of a conformally covariant <O T T> structure is
again conformally covariant (Delta_T = 3 = J_T + d - 2), and so is S(1,3) - S(3,1).  The conformal group acts on triples
of mutually spacelike points with a 1-dimensional stabiliser, so the orbit of one collinear triple is open in the
9-dimensional configuration space; both objects are real-analytic, hence they vanish identically iff they vanish at
ONE configuration for ALL polarisations.  We take x4 = 0, x3 = (0,0,1), x1 = (0,0,-2) (x1 kept symbolic until after
differentiation), z1 arbitrary (off the cone; the Todorov operator is interior), z3 = zeta(a3,1), z4 = zeta(a4,1)
(generic null polarisations up to scale; S is polynomial in z3 and, after removing V4^J, rational in z4).
The prefactor V4^J X13^{-e13} X14^{-e14} X34^{-e34} with SYMBOLIC exponents is differentiated analytically.
usage: python3 ap2_conservation.py [name ...]   (pp pm odd oddJ)"""
import pickle, sys, time
from common import *
import structures as st

x1 = list(sp.symbols('x1_0:3')); z1 = list(sp.symbols('z1_0:3'))
a1, a3, a4 = sp.symbols('a1 a3 a4')
ETA = (-1, 1, 1)
def zeta(a, b):
    return [a*b, (a*a + b*b)/2, sp.I*(a*a - b*b)/2]
X3 = [sp.Integer(0), sp.Integer(0), sp.Integer(1)]; X4 = [sp.Integer(0)]*3
X1v = [sp.Integer(0), sp.Integer(0), sp.Integer(-2)]
Z3, Z4 = zeta(a3, sp.Integer(1)), zeta(a4, sp.Integer(1))

def R_of(mr, B):
    m1, m3, k, n13, n14, n34, eps = mr
    r = B['V1']**m1*B['V3']**m3*B['V4']**(-k)*B['H13']**n13*B['H14']**n14*B['H34']**n34
    if eps:
        r = r*B['eps' + eps]*B['X' + eps]
    return r

def check_conservation(vec, parity):
    half = sp.Rational(1, 2) if parity == 'odd' else 0
    e13, e14 = 5 - S + half, S + half
    B = st.blocks(x1, z1, X3, Z3, X4, Z4)
    sub = dict(zip(x1, X1v))
    dlog = [(J*sp.diff(B['V4'], x1[m])/B['V4'] - e13*sp.diff(B['X13'], x1[m])/B['X13']
             - e14*sp.diff(B['X14'], x1[m])/B['X14']).subs(sub) for m in range(3)]
    tot = [0, 0, 0]
    for mr, c in vec.items():
        if c == 0:
            continue
        R = R_of(mr, B)
        dR = [sp.diff(R, z1[m]) for m in range(3)]
        box = sum(ETA[m]*sp.diff(dR[m], z1[m]) for m in range(3))
        DR = [sp.Rational(3, 2)*ETA[m]*dR[m] - sp.Rational(1, 2)*z1[m]*box for m in range(3)]
        div = sum(sp.diff(DR[m], x1[m]).subs(sub) + dlog[m]*DR[m].subs(sub) for m in range(3))
        for nu in range(3):     # div is linear in z1
            tot[nu] += c*sp.together(sp.diff(div, z1[nu]).subs({z: 0 for z in z1}))
    return [sp.expand(sp.numer(sp.together(t))) == 0 for t in tot]

def check_symmetry(vec, parity, Jpar):
    Z1 = zeta(a1, sp.Integer(1))
    B = st.blocks(X1v, Z1, X3, Z3, X4, Z4)
    Bs = st.blocks(X3, Z3, X1v, Z1, X4, Z4)
    assert sp.simplify(Bs['V4'] + B['V4']) == 0 and sp.simplify(Bs['X14'] - B['X34']) == 0 \
        and sp.simplify(Bs['X34'] - B['X14']) == 0 and sp.simplify(Bs['X13'] - B['X13']) == 0
    sgn = (-1)**Jpar          # prefactor symmetric up to V4^J -> (-V4)^J
    tot = sum(c*(R_of(mr, B) - sgn*R_of(mr, Bs)) for mr, c in vec.items() if c != 0)
    return sp.expand(sp.numer(sp.together(tot))) == 0

if __name__ == "__main__":
    structs = pickle.load(open(os.path.join(CACHE, 'symbolic_structures.pkl'), 'rb'))
    names = [] if sys.argv[1:2] == ['controls'] else (sys.argv[1:] or ['pp', 'pm', 'odd', 'oddJ'])
    for name in names:
        par = 'even' if name in ('pp', 'pm') else 'odd'
        t0 = time.time()
        cons = check_conservation(structs[name], par)
        t1 = time.time()
        sym = check_symmetry(structs[name], par, 1 if name == 'oddJ' else 0)
        print(f"{name:5s}: conservation (3 components) identically zero: {cons} [{t1-t0:.0f}s];  "
              f"Bose symmetry: {sym} [{time.time()-t1:.0f}s]", flush=True)

def controls():
    """negative controls: the test must detect (i) a perturbed coefficient, (ii) the pp structure with J -> J+2 in
    its coefficients (a wrong J-dependence), (iii) a symmetric non-conserved monomial combination."""
    structs = pickle.load(open(os.path.join(CACHE, 'symbolic_structures.pkl'), 'rb'))
    v = dict(structs['pp']); k0 = (1, 1, 0, 1, 0, 0, '')
    v1 = dict(v); v1[k0] = v1[k0] + 1
    v2 = {k: c.subs(J, J + 2) for k, c in v.items()}
    v3 = {(2, 2, 0, 0, 0, 0, ''): sp.Integer(1)}
    for nm, vv in (('perturbed', v1), ('J->J+2', v2), ('V1^2V3^2V4^J alone', v3)):
        print(f"control {nm}: conservation components zero? {check_conservation(vv, 'even')}", flush=True)

if __name__ == "__main__" and sys.argv[1:2] == ["controls"]:
    controls()
