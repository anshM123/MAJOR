"""
Embedding-space algebra (sympy) for spinning three- and four-point structures built with the differential
operators of Costa-Penedones-Poland-Rychkov (arXiv:1109.6321, sect. 3.3) acting on scalar structures.

Conventions.  P = (P^+, P^-, P^mu) = (1, x^2, x^mu),  Z = (0, 2 x.z, z^mu),  P.P' = -(x-x')^2/2.
Signature (-,+,+).  All symbols  s_ij := 2 P_i.P_j = -x_ij^2  (positive when x_ij is timelike, as in the frame
used for dispersive functionals),  z_ij := Z_i.P_j = z_i.(x_j - x_i),  z1z3 := Z_1.Z_3 = z_1.z_3.

Operators (spinning operator at point 1, scalar at point 2):
   D11 = (P1.P2)(Z1.d/dP2) - (Z1.P2)(P1.d/dP2)                              [Delta_1 -> Delta_1 - 1, spin_1 + 1]
   D12 = (P1.P2)(Z1.d/dP1) - (Z1.P2)(P1.d/dP1) + (Z1.P2)(Z1.d/dZ1)          [Delta_2 -> Delta_2 - 1, spin_1 + 1]
(the terms of 1109.6321 with d/dZ2 are absent because operator 2 is a scalar).  Both are tangent to P^2 = 0, Z.P = 0
and preserve transversality (checked in t2_structures.py).  D33, D34: the same with (1,2) -> (3,4).
"""
import sympy as sp, itertools, pickle, os

_here = os.path.dirname(os.path.abspath(__file__))

class Engine:
    """Functions of dot products of a fixed list of embedding vectors; derivations V.d/dW."""
    def __init__(self, dots):
        # dots: dict {(A,B) sorted tuple: symbol or 0}
        self.dots = dots
        self.syms = {v: k for k, v in dots.items() if v != 0}
        self.logK = None          # log of the prefactor (sympy expr in the symbols)
        self.G = None             # dict (i,j) -> symbol for d_U^i d_V^j g
        self.U = self.V = None

    def dot(self, A, B):
        return self.dots.get(tuple(sorted((A, B))), 0)

    def dsym(self, q, V, W):
        A, B = self.syms[q]
        r = 0
        if A == W: r += self.dot(V, B)
        if B == W: r += self.dot(A, V)
        return r

    def dexpr(self, e, V, W):
        r = 0
        for q in e.free_symbols:
            if q in self.syms:
                d = self.dsym(q, V, W)
                if d != 0:
                    r += sp.diff(e, q)*d
        return r

    def dirK(self, Phi, V, W):
        """(1/K) V.d/dW (K Phi)"""
        r = self.dexpr(Phi, V, W) + self.dexpr(self.logK, V, W)*Phi
        if self.G is not None:
            dU = self.dexpr(self.U, V, W); dV = self.dexpr(self.V, V, W)
            for (i, j), g in self.G.items():
                c = sp.diff(Phi, g)
                if c != 0:
                    if (i+1, j) not in self.G or (i, j+1) not in self.G:
                        raise ValueError("jet order exceeded")
                    r += c*(self.G[(i+1, j)]*dU + self.G[(i, j+1)]*dV)
        return r

    def Dop(self, Phi, kind, pt=(1, 2)):
        """kind 'a': D_{pp} (shifts the dimension at the spinning point), 'b': D_{pq} (shifts the scalar's)."""
        p, q = 'P%d' % pt[0], 'P%d' % pt[1]
        Zp = 'Z%d' % pt[0]
        PP, ZP = self.dot(p, q), self.dot(Zp, q)
        if kind == 'a':
            return sp.expand(PP*self.dirK(Phi, Zp, q) - ZP*self.dirK(Phi, p, q))
        # Z.d/dZ acts only on Phi (the prefactor K contains no Z of this point)
        return sp.expand(PP*self.dirK(Phi, Zp, p) - ZP*self.dirK(Phi, p, p) + ZP*self.dexpr(Phi, Zp, Zp))


def engine3():
    """three-point: operators 1 (spinning, Z1), 2 (scalar), 5 (spin J, polarisation Z5)."""
    s12, s15, s25 = sp.symbols('s12 s15 s25', positive=True)
    z12, z15, z51, z52, z1z5 = sp.symbols('z12 z15 z51 z52 z1z5')
    d = {('P1', 'P2'): s12/2, ('P1', 'P5'): s15/2, ('P2', 'P5'): s25/2,
         ('P2', 'Z1'): z12, ('P5', 'Z1'): z15, ('P1', 'Z5'): z51, ('P2', 'Z5'): z52, ('Z1', 'Z5'): z1z5}
    # Engine.syms needs symbol -> pair; s_ij/2 entries: register the bare symbols
    e = Engine(d)
    e.syms = {s12: ('P1', 'P2'), s15: ('P1', 'P5'), s25: ('P2', 'P5'), z12: ('P2', 'Z1'), z15: ('P5', 'Z1'),
              z51: ('P1', 'Z5'), z52: ('P2', 'Z5'), z1z5: ('Z1', 'Z5')}
    e._half = {s12, s15, s25}
    return e, dict(s12=s12, s15=s15, s25=s25, z12=z12, z15=z15, z51=z51, z52=z52, z1z5=z1z5)

def _patch_dsym(e):
    """d(s_ij) = 2 d(P_i.P_j)."""
    old = e.dsym
    def dsym(q, V, W):
        r = old(q, V, W)
        return 2*r if q in e._half else r
    e.dsym = dsym

Delta, J, D1, D2, D3, D4 = sp.symbols('Delta J Delta1 Delta2 Delta3 Delta4')

def seed3():
    """scalar-scalar-spin J structure  S3 = V5^J / (s12^al s15^be s25^ga),  V5 = (Z5.P1)(P2.P5) - (Z5.P2)(P1.P5).
    In position space this is  [ (n.x15) x25^2 - (n.x25) x15^2 ]^J (-1/2)^J ... / (|x12|^{..} ...)   (see t2_structures.py)."""
    e, sy = engine3(); _patch_dsym(e)
    V5 = (sy['z51']*sy['s25'] - sy['z52']*sy['s15'])/2
    al = (D1 + D2 - Delta + J)/2; be = (Delta + J + D1 - D2)/2; ga = (Delta + J - D1 + D2)/2
    e.logK = J*sp.log(V5) - al*sp.log(sy['s12']) - be*sp.log(sy['s15']) - ga*sp.log(sy['s25'])
    return e, sy, V5

def word3(word):
    """Q = (D-word S3)/S3 for the seed with dimensions (Delta1 + #a, Delta2 + #b); word: string of 'a','b' applied right to left."""
    e, sy, V5 = seed3()
    na, nb = word.count('a'), word.count('b')
    e.logK = e.logK.subs({D1: D1 + na, D2: D2 + nb}, simultaneous=True)
    Phi = sp.Integer(1)
    for kind in reversed(word):
        Phi = e.Dop(Phi, kind, (1, 2))
    return sp.factor(sp.simplify(Phi)), sy, V5

def engine4():
    names = ['s12', 's13', 's14', 's23', 's24', 's34', 'z12', 'z13', 'z14', 'z31', 'z32', 'z34', 'z1z3']
    S = {n: sp.Symbol(n, positive=n.startswith('s')) for n in names}
    d = {}
    pairs = {'s12': ('P1', 'P2'), 's13': ('P1', 'P3'), 's14': ('P1', 'P4'), 's23': ('P2', 'P3'), 's24': ('P2', 'P4'), 's34': ('P3', 'P4'),
             'z12': ('P2', 'Z1'), 'z13': ('P3', 'Z1'), 'z14': ('P4', 'Z1'), 'z31': ('P1', 'Z3'), 'z32': ('P2', 'Z3'), 'z34': ('P4', 'Z3'),
             'z1z3': ('Z1', 'Z3')}
    for n, pr in pairs.items():
        d[tuple(sorted(pr))] = S[n]/2 if n.startswith('s') else S[n]
    e = Engine(d)
    e.syms = {S[n]: tuple(sorted(pr)) for n, pr in pairs.items()}
    e._half = {S[n] for n in names if n.startswith('s')}
    _patch_dsym(e)
    return e, S

JETMAX = 4
def Gsyms(order=JETMAX):
    return {(i, j): sp.Symbol('G_%d_%d' % (i, j)) for i in range(order+1) for j in range(order+1-i)}

def word4(wordL, wordR):
    """Phi with  D-words acting on the scalar partial wave written as
          K g^{a,b}(u,v) = (-1)^J Kt gt(P, Q),     gt := g^{-a,b}_{Delta,J} at (u/v, 1/v)   [g^{a,b}(u,v) = (-1)^J v^{-b} g^{-a,b}(u/v,1/v)],
          Kt = K v^{-b} = (s24/s14)^{(d1-d2)/2} (s24/s23)^{(d3-d4)/2} / (s12^{(d1+d2)/2} s34^{(d3+d4)/2}),
          P = (log u + log v)/2,  Q = (log v - log u)/2,   u = s12 s34/(s13 s24), v = s14 s23/(s13 s24).
    The result of the D-words is (-1)^J Kt * Phi(s, z; G_ij),  G_ij = d_P^i d_Q^j gt.  Seed dimensions:
    (d1,d2,d3,d4) = (D1 + #a_L, D2 + #b_L, D3 + #a_R, D4 + #b_R).
    (The variables P, Q and the reflected block are used because they are well conditioned in the Regge region:
    d_P is O(1) while d_Q is large, and gt is a sum of positive terms for rho, rhobar < 0.)"""
    e, S = engine4()
    d1, d2 = D1 + wordL.count('a'), D2 + wordL.count('b')
    d3, d4 = D3 + wordR.count('a'), D4 + wordR.count('b')
    lg = {k: sp.log(v) for k, v in S.items() if k.startswith('s')}
    e.logK = ((d1-d2)/2*(lg['s24'] - lg['s14']) + (d3-d4)/2*(lg['s24'] - lg['s23'])
              - (d1+d2)/2*lg['s12'] - (d3+d4)/2*lg['s34'])
    U = lg['s12'] + lg['s34'] - lg['s13'] - lg['s24']
    V = lg['s14'] + lg['s23'] - lg['s13'] - lg['s24']
    e.U = (U + V)/2          # P
    e.V = (V - U)/2          # Q
    e.G = Gsyms()
    Phi = e.G[(0, 0)]
    for kind in reversed(wordL):
        Phi = e.Dop(Phi, kind, (1, 2))
    for kind in reversed(wordR):
        Phi = e.Dop(Phi, kind, (3, 4))
    return Phi, S, e.G

def cached(name, fun):
    fn = os.path.join(_here, 'cache_' + name + '.pkl')
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    r = fun()
    pickle.dump(r, open(fn, 'wb'))
    return r
