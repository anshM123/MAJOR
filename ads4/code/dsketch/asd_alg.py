"""
asd_alg: truncated two-parameter perturbation algebra on Euclidean Poincare AdS4 (DSKETCH_UP.md, Target 1(a)).

A field is a dict {(p, q): expr(z)} with p, q in {0, 1}; the term of order a^p b^q carries the plane-wave factor
exp(i (p k3 + q k4) . x), which is NOT written: d/dx_j acts as multiplication by i (p k3 + q k4)_j and d/dz as d/dz.
Products are truncated at a^2 = b^2 = 0, so (1,1) is the a*b ("second-order, cross") term.
Coordinates X = (x1, x2, x3, z); the 4D orientation is eps_{x1 x2 x3 z} = +1 (same as c1val/v2_ads4_sd_modes.py).
"""
import itertools
import sympy as sp

z = sp.symbols('z', positive=True)
N = 4
ORDERS = [(0, 0), (1, 0), (0, 1), (1, 1)]


class Alg:
    def __init__(self, k3, k4):
        self.k3 = [sp.nsimplify(c) for c in k3]
        self.k4 = [sp.nsimplify(c) for c in k4]

    def mom(self, o):
        return [o[0]*self.k3[j] + o[1]*self.k4[j] for j in range(3)]

    # --- series operations -------------------------------------------------------------------------------------------
    @staticmethod
    def add(*fs):
        out = {}
        for f in fs:
            for o, e in f.items():
                out[o] = out.get(o, 0) + e
        return out

    @staticmethod
    def scale(c, f):
        return {o: c*e for o, e in f.items()}

    @staticmethod
    def mul(f, g):
        out = {}
        for o1, e1 in f.items():
            if e1 == 0:
                continue
            for o2, e2 in g.items():
                if e2 == 0:
                    continue
                o = (o1[0] + o2[0], o1[1] + o2[1])
                if o[0] > 1 or o[1] > 1:
                    continue
                out[o] = out.get(o, 0) + e1*e2
        return out

    def d(self, f, mu):
        out = {}
        for o, e in f.items():
            if mu == 3:
                out[o] = sp.diff(e, z)
            else:
                out[o] = sp.I*self.mom(o)[mu]*e
        return out

    @staticmethod
    def get(f, o):
        return f.get(o, 0)

    @staticmethod
    def simp(f):
        return {o: sp.simplify(e) for o, e in f.items()}


def christoffel(A, g, gi):
    dg = [[[A.d(g[b][c], e) for e in range(N)] for c in range(N)] for b in range(N)]   # dg[b][c][e] = d_e g_bc
    G = [[[None]*N for _ in range(N)] for _ in range(N)]
    for a, b, c in itertools.product(range(N), repeat=3):
        acc = {}
        for e in range(N):
            t = A.add(dg[e][b][c], dg[e][c][b], A.scale(-1, dg[b][c][e]))
            acc = A.add(acc, A.mul(gi[a][e], t))
        G[a][b][c] = A.scale(sp.Rational(1, 2), acc)
    return G


def riemann_up(A, G):
    """R^a_{bcd} = d_c G^a_{db} - d_d G^a_{cb} + G^a_{ce} G^e_{db} - G^a_{de} G^e_{cb}"""
    R = {}
    for a, b, c, dd in itertools.product(range(N), repeat=4):
        t = A.add(A.d(G[a][dd][b], c), A.scale(-1, A.d(G[a][c][b], dd)))
        for e in range(N):
            t = A.add(t, A.mul(G[a][c][e], G[e][dd][b]), A.scale(-1, A.mul(G[a][dd][e], G[e][c][b])))
        R[(a, b, c, dd)] = t
    return R


def ricci(A, Rup):
    return [[A.add(*[Rup[(a, b, a, dd)] for a in range(N)]) for dd in range(N)] for b in range(N)]


def inverse_conformal(A, H):
    """inverse of (delta + H) to the truncated order: delta - H + H.H (only the (1,1) part of H.H survives)."""
    gi = [[None]*N for _ in range(N)]
    for i, j in itertools.product(range(N), repeat=2):
        t = {(0, 0): sp.Integer(1 if i == j else 0)}
        t = A.add(t, A.scale(-1, H[i][j]))
        for e in range(N):
            t = A.add(t, A.mul(H[i][e], H[e][j]))
        gi[i][j] = t
    return gi


def helicity_basis(kvec):
    """unit khat, and right-handed (n1, n2, khat); e^{+-} = n1 +- i n2 (the labels of c1val/v2: k along x3 -> n1 = x1, n2 = x2)."""
    k = sp.Matrix([sp.nsimplify(c) for c in kvec])
    kn = sp.sqrt(k.dot(k))
    kh = k/kn
    # pick n1 perpendicular to k: try x2, then x1
    for trial in (sp.Matrix([0, 1, 0]), sp.Matrix([1, 0, 0]), sp.Matrix([0, 0, 1])):
        n1 = trial - trial.dot(kh)*kh
        if n1.dot(n1) != 0:
            break
    n1 = n1/sp.sqrt(n1.dot(n1))
    n2 = kh.cross(n1)
    assert sp.simplify(n1.cross(n2) - kh) == sp.zeros(3, 1)
    return kn, kh, n1, n2


def mode(kvec, hel, regular=True):
    """TT radial-gauge graviton mode gamma_ij(z) (conformally rescaled: h = gamma/z^2) without the plane-wave factor."""
    kn, kh, n1, n2 = helicity_basis(kvec)
    e = n1 + hel*sp.I*n2
    phi = (1 + kn*z)*sp.exp(-kn*z) if regular else (1 - kn*z)*sp.exp(kn*z)
    gam = sp.zeros(4, 4)
    for i in range(3):
        for j in range(3):
            gam[i, j] = sp.simplify(e[i]*e[j]*phi)
    return gam
