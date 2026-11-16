"""
Embedding-space algebra for three- and four-point structures with TWO spinning operators per vertex, built with the
differential operators of Costa-Penedones-Poland-Rychkov (arXiv:1109.6321, sect. 3) acting on scalar seeds.

Conventions as ads4/code/e2e/embed.py:  P = (1, x^2, x),  Z = (0, 2 x.z, z),  P.P' = -(x-x')^2/2, signature (-,+,+);
symbols  s_ij := 2 P_i.P_j,  z_ij := Z_i.P_j,  zizj := Z_i.Z_j.

Operators at a vertex (i, j) (both points may be spinning):
   D_ii = (P_i.P_j)(Z_i.d/dP_j) - (Z_i.P_j)(P_i.d/dP_j) - (Z_i.Z_j)(P_i.d/dZ_j) + (P_i.Z_j)(Z_i.d/dZ_j)
          [spin at i +1, Delta_i -> Delta_i - 1]
   D_ij = (P_i.P_j)(Z_i.d/dP_i) - (Z_i.P_j)(P_i.d/dP_i) + (Z_i.P_j)(Z_i.d/dZ_i)
          [spin at i +1, Delta_j -> Delta_j - 1]
   H_ij = (Z_i.Z_j)(P_i.P_j) - (Z_i.P_j)(Z_j.P_i)            [multiplication; spins at i and j +1, both Delta's -1 ... see seeds]
(for a scalar at j the Z_j terms vanish and D_ii, D_ij are the operators 'a', 'b' of e2e/embed.py).  The operators annihilate
the constraints P_k^2, Z_k.P_k, Z_k^2 (checked in t0_structures.py), so the dot products can be treated as independent symbols.

Word notation (string, operators applied RIGHT TO LEFT):  'A' = D_11, 'B' = D_12, 'C' = D_22, 'E' = D_21 (vertex (1,2));
for the vertex (3,4) the same letters denote D_33, D_34, D_44, D_43;  'H' = H_12 (resp. H_34).
Seed shifts: A, E, H raise the seed dimension at the first point of the vertex by one; B, C, H raise the second one.
"""
import sympy as sp, os, pickle, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
from embed import Engine            # noqa: E402  (generic derivation engine on dot products)

Delta, J, D1, D2, D3, D4 = sp.symbols('Delta J Delta1 Delta2 Delta3 Delta4')
CACHE = os.path.join(HERE, 'cache')

def shifts(word):
    """(shift of the first point's dimension, shift of the second point's dimension) of the seed for this word"""
    return (word.count('A') + word.count('E') + word.count('H'), word.count('B') + word.count('C') + word.count('H'))

def _mk_engine(names_pairs):
    S = {n: sp.Symbol(n, positive=n.startswith('s')) for n in names_pairs}
    d = {}
    for n, pr in names_pairs.items():
        d[tuple(sorted(pr))] = S[n]/2 if n.startswith('s') else S[n]
    e = Engine(d)
    e.syms = {S[n]: tuple(sorted(pr)) for n, pr in names_pairs.items()}
    half = {S[n] for n in names_pairs if n.startswith('s')}
    old = e.dsym
    def dsym(q, V, W):
        r = old(q, V, W)
        return 2*r if q in half else r
    e.dsym = dsym
    return e, S

def op(e, Phi, kind, i, j):
    """apply one operator at vertex (i, j) to K*Phi, return the new Phi (prefactor K = exp(e.logK) unchanged)."""
    Pi, Pj, Zi, Zj = 'P%d' % i, 'P%d' % j, 'Z%d' % i, 'Z%d' % j
    PP = e.dot(Pi, Pj)
    if kind == 'H':
        return sp.expand(Phi*(e.dot(Zi, Zj)*PP - e.dot(Zi, Pj)*e.dot(Zj, Pi)))
    if kind in ('A', 'C'):
        if kind == 'C':
            Pi, Pj, Zi, Zj = Pj, Pi, Zj, Zi
        r = PP*e.dirK(Phi, Zi, Pj) - e.dot(Zi, Pj)*e.dirK(Phi, Pi, Pj)
        r += -e.dot(Zi, Zj)*e.dexpr(Phi, Pi, Zj) + e.dot(Pi, Zj)*e.dexpr(Phi, Zi, Zj)
        return sp.expand(r)
    if kind in ('B', 'E'):
        if kind == 'E':
            Pi, Pj, Zi, Zj = Pj, Pi, Zj, Zi
        r = PP*e.dirK(Phi, Zi, Pi) - e.dot(Zi, Pj)*e.dirK(Phi, Pi, Pi) + e.dot(Zi, Pj)*e.dexpr(Phi, Zi, Zi)
        return sp.expand(r)
    raise ValueError(kind)

# ------------------------------------------------------------------ three points: O_1 (Z1), O_2 (Z2), exchanged O_5 (Z5)
PAIRS3 = {'s12': ('P1', 'P2'), 's15': ('P1', 'P5'), 's25': ('P2', 'P5'),
          'z12': ('P2', 'Z1'), 'z15': ('P5', 'Z1'), 'z21': ('P1', 'Z2'), 'z25': ('P5', 'Z2'), 'z51': ('P1', 'Z5'), 'z52': ('P2', 'Z5'),
          'z1z2': ('Z1', 'Z2'), 'z1z5': ('Z1', 'Z5'), 'z2z5': ('Z2', 'Z5')}

def word3(word):
    """Q_w := (word applied to S3(D1 + sh1, D2 + sh2)) / S3(D1 + sh1, D2 + sh2) as a polynomial/rational expression,
    S3(d1, d2) = v5^J / (s12^al s15^be s25^ga),  v5 = (Z5.P1)(P2.P5) - (Z5.P2)(P1.P5),
    al = (d1+d2-Delta+J)/2, be = (Delta+J+d1-d2)/2, ga = (Delta+J-d1+d2)/2.   Returns (Q, symbols, v5, (d1, d2))."""
    e, S = _mk_engine(PAIRS3)
    v5 = (S['z51']*S['s25'] - S['z52']*S['s15'])/2
    sh1, sh2 = shifts(word)
    d1, d2 = D1 + sh1, D2 + sh2
    al = (d1 + d2 - Delta + J)/2; be = (Delta + J + d1 - d2)/2; ga = (Delta + J - d1 + d2)/2
    e.logK = J*sp.log(v5) - al*sp.log(S['s12']) - be*sp.log(S['s15']) - ga*sp.log(S['s25'])
    Phi = sp.Integer(1)
    for kind in reversed(word):
        Phi = op(e, Phi, kind, 1, 2)
    return sp.factor(sp.together(Phi)), S, v5, (d1, d2)

# ------------------------------------------------------------------ four points
SN = ['s12', 's13', 's14', 's23', 's24', 's34']
ZN = ['z12', 'z13', 'z14', 'z21', 'z23', 'z24', 'z31', 'z32', 'z34', 'z41', 'z42', 'z43', 'z1z2', 'z1z3', 'z1z4', 'z2z3', 'z2z4', 'z3z4']

def pairs4():
    d = {}
    for n in SN:
        d[n] = ('P' + n[1], 'P' + n[2])
    for n in ZN:
        if n.count('z') == 2:
            d[n] = ('Z' + n[1], 'Z' + n[3])
        else:
            d[n] = ('Z' + n[1], 'P' + n[2])
    return d

def Gsyms(order):
    return {(i, j): sp.Symbol('G_%d_%d' % (i, j)) for i in range(order + 1) for j in range(order + 1 - i)}

def word4(wordL, wordR, jet):
    """Phi with  D-words acting on the scalar partial wave  K g^{a,b}(u,v) = (-1)^J Kt gt(P, Q)  (conventions of e2e/embed.word4):
       Kt = (s24/s14)^{(d1-d2)/2} (s24/s23)^{(d3-d4)/2} / (s12^{(d1+d2)/2} s34^{(d3+d4)/2}),
       gt = g^{-a,b}(u/v, 1/v),  P = (log u + log v)/2,  Q = (log v - log u)/2,  G_ij = d_P^i d_Q^j gt.
    Seed dimensions (d1,d2,d3,d4) = (D1 + shL1, D2 + shL2, D3 + shR1, D4 + shR2).  Result: (-1)^J Kt Phi."""
    e, S = _mk_engine(pairs4())
    s1, s2 = shifts(wordL); s3, s4 = shifts(wordR)
    d1, d2, d3, d4 = D1 + s1, D2 + s2, D3 + s3, D4 + s4
    lg = {k: sp.log(v) for k, v in S.items() if k.startswith('s')}
    e.logK = ((d1 - d2)/2*(lg['s24'] - lg['s14']) + (d3 - d4)/2*(lg['s24'] - lg['s23'])
              - (d1 + d2)/2*lg['s12'] - (d3 + d4)/2*lg['s34'])
    U = lg['s12'] + lg['s34'] - lg['s13'] - lg['s24']
    V = lg['s14'] + lg['s23'] - lg['s13'] - lg['s24']
    e.U = (U + V)/2
    e.V = (V - U)/2
    e.G = Gsyms(jet)
    Phi = e.G[(0, 0)]
    for kind in reversed(wordL):
        Phi = op(e, Phi, kind, 1, 2)
    for kind in reversed(wordR):
        Phi = op(e, Phi, kind, 3, 4)
    return Phi, S, e.G

def cached(name, fun):
    os.makedirs(CACHE, exist_ok=True)
    fn = os.path.join(CACHE, name + '.pkl')
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    r = fun()
    pickle.dump(r, open(fn, 'wb'))
    return r
