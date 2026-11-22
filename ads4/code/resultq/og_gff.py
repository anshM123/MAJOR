"""
og_gff (RESULT_Q_OG_CHECK.md §3.5): the same swapped-form Lemma-A test as og_free, for a one-parameter family of spinning correlators with
GENERIC (non-half-integer) exponents: O = [phi phi]_{0, J=2}, the leading-twist spin-2 double-trace primary of a generalised free scalar
phi of dimension Dp (Delta_O = 2 Dp + 2, J = 2), connected (O(1/N)) part of <O O O O>.  Dp = 1/2 is the free-scalar T of og_free.

* links: (sigma X^2)^{-Dp} with the same branch / ^+ rules as og_free (mixed links positive on the real domain);
* vertex O(x, z) = a0 [phi (z.d)^2 phi] + a1 (z.d phi)^2, with a0/a1 fixed by requiring a conformal 2-pt function (checked);
* light-transform weight t^(Delta_O + J - 2) = t^(2 Dp + 2).
"""
import itertools
import numpy as np
import sympy as sp
import og_free as F

dot, zvec, E0 = F.dot, F.zvec, F.E0


def _derive(MMAX=2):
    r, Y, A, B, ab, sg, p = sp.symbols('r Y A B ab sg p')

    def dY(e):
        return sp.diff(e, Y) + sp.diff(e, r)*(-p*r/Y)

    def da(e):
        return dY(e)*2*sg*A + sp.diff(e, B)*ab

    def db(e):
        return dY(e)*2*sg*B + sp.diff(e, A)*ab
    out = {}
    for n in range(MMAX + 1):
        e = r
        for _ in range(n):
            e = -db(e)
        for m in range(MMAX + 1):
            em = e
            for _ in range(m):
                em = da(em)
            out[(m, n)] = sp.lambdify((r, Y, A, B, ab, sg, p), sp.simplify(em), 'numpy')
    return out


DFUN = _derive()


def twopt_ratio(Dp, a0, a1, X, z, w):
    """<O(x,z) O(0,w)> / [((I z).w)^2 / (X^2)^{2 Dp + 2}] for a real spacelike X."""
    X2 = dot(X, X); r = X2**-Dp
    Dm = lambda m, n: DFUN[(m, n)](r, X2, dot(z, X), dot(w, X), dot(z, w), 1.0, Dp)
    V = [(a1, (1, 1)), (a0, (0, 2))]
    tot = sum(c1*c2*(Dm(a, c)*Dm(b, d) + Dm(a, d)*Dm(b, c)) for c1, (a, b) in V for c2, (c, d) in V)
    Iz = z - 2*X*dot(X, z)/X2
    return tot/(dot(Iz, w)**2/X2**(2*Dp + 2))


def primary_coef(Dp, seed=3):
    """a0/a1 of the spin-2 leading-twist double-trace primary: the conformality condition on the 2-pt function is a perfect square in
    a0 with the double root a0 = -Dp/(Dp + 1) (found numerically, identified, and checked here; Dp = 1/2: -1/3 = improved free T).
    Returns (a0, max relative spread of the 2-pt ratio over 6 random spacelike configurations)."""
    a0 = -Dp/(Dp + 1.0)
    rng = np.random.default_rng(seed)
    rs = []
    for _ in range(6):
        th = rng.uniform(0, 2*np.pi); X = rng.uniform(0.7, 1.5)*np.array([rng.uniform(-0.5, 0.5), np.cos(th), np.sin(th)])
        z = zvec(rng.uniform(0, 6))*rng.uniform(0.5, 2); w = zvec(rng.uniform(0, 6))*rng.uniform(0.5, 2)
        rs.append(twopt_ratio(Dp, a0, 1.0, X, z, w))
    rs = np.array(rs)
    return a0, float(np.max(np.abs(rs - rs[0]))/abs(rs[0]))


def links(t1, t3, z1, z3, pol2, pol4, y13, Dp):
    t1 = np.asarray(t1, complex); t3 = np.asarray(t3, complex)
    sh = np.broadcast(t1, t3).shape
    t1 = np.broadcast_to(t1, sh); t3 = np.broadcast_to(t3, sh)
    X = {1: -t1[..., None]*z1, 2: np.broadcast_to(-E0, sh + (3,)).astype(complex),
         3: t3[..., None]*z3, 4: np.broadcast_to(E0, sh + (3,)).astype(complex)}
    P = {1: z1, 2: pol2, 3: z3, 4: pol4}
    out = {}
    for (i, j) in F.PAIRS:
        sg = -1.0 if (i, j) in F.MIXED else 1.0
        if (i, j) == (1, 3):
            Y = 2*y13*t1*t3
            r = (2*y13)**-Dp*t1**-Dp*t3**-Dp
            A = t3*y13; B = t1*y13; ab = -y13
        else:
            Xij = X[i] - X[j]
            if (i, j) == (1, 2):
                Y = 2*t1 - 1 + 0j
            elif (i, j) == (3, 4):
                Y = 2*t3 - 1 + 0j
            else:
                Y = sg*dot(Xij, Xij)
            r = Y**-Dp
            A = dot(P[i], Xij); B = dot(P[j], Xij); ab = float(dot(P[i], P[j]))
        out[(i, j)] = {mn: DFUN[mn](r, Y, A, B, ab, sg, Dp) for mn in DFUN}
    return out


def Eplus(k, phi1, phi3, pol2, pol4, Dp, a0, n=80, rad=0.5, classes=('1234', '1243', '1324')):
    z1 = zvec(phi1); z3 = zvec(phi3); y13 = 2*np.sin((phi1 - phi3)/2)**2
    t, w = F.contour(n, rad)
    T1 = t[:, None]; T3 = t[None, :]
    L = links(T1, T3, z1, z3, pol2, pol4, y13, Dp)
    keep = list(F.VERTEX)
    F.VERTEX[:] = [(1.0, (1, 1)), (a0, (0, 2))]
    try:
        G = F.correlator(L, classes)
    finally:
        F.VERTEX[:] = keep
    T1b, T3b = np.broadcast_arrays(T1, T3)
    h0 = T1b**(2*Dp + 2)*T3b**(2*Dp + 2)*F.fk(k, T1b, T3b, z1, z3, y13=y13)
    W = w[:, None]*w[None, :]
    return {cl: complex(np.sum(W*h0*G[cl])) for cl in classes}
