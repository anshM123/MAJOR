"""
Exact mean-field-theory OPE coefficients of the double-twist primaries [O phi]_{N,J} from the two-particle Fock space
(no conformal blocks used).  d = 3, Euclidean.

  O : conserved, spin l (l = 1: Delta_O = 2;  l = 2: Delta_O = 3)   [l = 0: a scalar of dimension D1, used as a test]
  phi : scalar of dimension Dphi (rational)

Method (derivation in SPINNING_MFT.md sect. 2):
 * one-particle states are represented by polynomial wavefunctions f(x, z) = <f| O(x, z) |0>; the inner product is the
   one whose reproducing kernel is  K(y, w; x, z) = <0| O(y, w)^dagger O(x, z) |0>,
       K_phi = D^{-Dphi},   K_O = (w.I(y)I(y/y^2 - x).z)^l / D^{Delta_O} = A(y,w;x,z)^l / D^{2l+1},   D = 1 - 2x.y + x^2 y^2
   (A is a polynomial, bi-degree <= (2,2)).  If K = sum_{alpha,beta} y^alpha M_{alpha beta} x^beta (polarisations included,
   reduced modulo z^2 = w^2 = 0), the norm of a wavefunction psi is psi^T M^+ conj(psi) (M^+ on the range of M).
 * a primary X of dimension Delta = Delta_O + Dphi + N and spin J has two-particle wavefunction
       q_X(x1, x2; z1, zeta) = <0| X(inf, zeta) O(x1, z1) phi(x2) |0>,   X(inf, zeta) := lim L^{2 Delta} X(L e, I(e) zeta),
   which is translation invariant (P annihilates <X|), i.e. a polynomial 3-pt structure T(x1 - x2), conserved in x1.
   Conversely every such polynomial structure is a primary.  <X,zbar|X,zbar> = N_X (zeta.zbar)^J.
 * with <O1 phi2 X> = lam * (structure S) for unit-normalised X:   lam^2 = (zeta.zbar)^J / <S(inf), S(inf)>_Fock.
All arithmetic is exact (python-flint rationals).
"""
import itertools, functools
import sympy as sp
from flint import fmpq, fmpq_mat

X = sp.symbols('x0:3'); Y = sp.symbols('y0:3'); Z = sp.symbols('z0:3'); W = sp.symbols('w0:3')
xv, yv, zv, wv = (sp.Matrix(v) for v in (X, Y, Z, W))
dot = lambda a, b: (a.T*b)[0]

def monos(a):
    return [(i, j, a-i-j) for i in range(a, -1, -1) for j in range(a-i, -1, -1)]

def pol_basis(l):
    """basis of polynomials of degree l in the polarisation, modulo z.z (Euclidean): exponent tuples"""
    if l == 0: return [(0, 0, 0)]
    if l == 1: return [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    if l == 2: return [(2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1)]
    raise ValueError

def reduce_null(expr, v, l):
    """reduce a polynomial homogeneous of degree l<=2 in v modulo v.v (replace v2^2 -> -v0^2 - v1^2)."""
    if l < 2: return sp.expand(expr)
    return sp.expand(sp.expand(expr).subs(v[2]**2, -v[0]**2 - v[1]**2))

@functools.lru_cache(None)
def A_poly():
    D = 1 - 2*dot(xv, yv) + dot(xv, xv)*dot(yv, yv)
    w_ = yv/dot(yv, yv) - xv
    I = lambda v: sp.eye(3) - 2*v*v.T/dot(v, v)
    return sp.expand(sp.cancel(sp.together((wv.T*I(yv)*I(w_)*zv)[0]*D)))

def bideg_parts(expr):
    """split a polynomial in (X,Y) into parts of bi-degree (a,a); returns dict a -> expr"""
    P = sp.Poly(expr, *X)
    out = {}
    for mon, c in P.terms():
        a = sum(mon)
        out[a] = out.get(a, 0) + c*sp.prod([X[i]**mon[i] for i in range(3)])
    return {a: sp.expand(v) for a, v in out.items()}

def Dinv_part(s, a):
    """bi-degree-(a,a) part of D^{-s}, D = 1 - t, t = 2x.y - x^2 y^2"""
    xy, xx, yy = dot(xv, yv), dot(xv, xv), dot(yv, yv)
    tot = 0
    for k in range(a + 1):
        j = a - k
        if j > k: continue
        tot += sp.rf(s, k)/sp.factorial(k)*sp.binomial(k, j)*(2*xy)**(k-j)*(-xx*yy)**j
    return sp.expand(tot)

def kernel_matrix(kind, a, l=0, s=None):
    """M^{(a)}: rows = (y-monomial, w-pol), cols = (x-monomial, z-pol).  kind 'scalar' (dimension s) or 'cons' (spin l)."""
    if kind == 'scalar':
        K = Dinv_part(s, a)
        pb = [(0, 0, 0)]
    else:
        Al = sp.expand(A_poly()**l)
        Aparts = bideg_parts(Al)
        K = 0
        for i, Ai in Aparts.items():
            if i <= a:
                K += sp.expand(Ai*Dinv_part(2*l + 1, a - i))
        K = reduce_null(reduce_null(K, Z, l), W, l)
        pb = pol_basis(l)
    ms = monos(a)
    rows = [(m, p) for m in ms for p in pb]
    idx = {r: i for i, r in enumerate(rows)}
    n = len(rows)
    M = fmpq_mat(n, n)
    P = sp.Poly(K, *(Y + W + X + Z))
    for mon, c in P.terms():
        r = idx[(mon[0:3], mon[3:6])]; cc = idx[(mon[6:9], mon[9:12])]
        cq = sp.Rational(c)
        M[r, cc] = fmpq(int(cq.p), int(cq.q))
    return M, rows

def _rref_pivots(M):
    R, rank = M.rref()
    piv = []
    r = 0
    for c in range(M.ncols()):
        if r < rank and R[r, c] != 0:
            piv.append(c); r += 1
    return piv

@functools.lru_cache(None)
def gram_inverse(kind, a, l=0, s=None):
    """returns (rows, pivot list P, inverse of M_PP)"""
    M, rows = kernel_matrix(kind, a, l, s)
    piv = _rref_pivots(M)
    n = len(piv)
    MP = fmpq_mat(n, n)
    for i, p in enumerate(piv):
        for j, q in enumerate(piv):
            MP[i, j] = M[p, q]
    return rows, piv, MP.inv(), M

def wavefunction_matrices(q, l, N):
    """q: sympy polynomial in X (= x1 - x2), Z (O polarisation), with complex numeric coefficients.
    returns dict a -> (Qr, Qi, rowsO, rowsPhi) where Q[(x1-mono, zpol), x2-mono] for deg_{x1} = a."""
    X1 = sp.symbols('u0:3'); X2 = sp.symbols('v0:3')
    qq = sp.expand(q.subs({X[i]: X1[i] - X2[i] for i in range(3)}, simultaneous=True))
    qq = reduce_null(qq, Z, l)
    P = sp.Poly(qq, *(X1 + X2 + Z))
    pb = pol_basis(l)
    data = {}
    for mon, c in P.terms():
        a = sum(mon[0:3])
        data.setdefault(a, {})[(mon[0:3], mon[6:9], mon[3:6])] = c
    return data

def fock_norm(q, l, N, Dphi, D1=None):
    """<q, q>_Fock for a two-particle wavefunction q(x1 - x2) of total degree N (exact)."""
    data = wavefunction_matrices(q, l, N)
    s = sp.Rational(Dphi)
    tot = sp.Integer(0)
    for a, entries in data.items():
        if l == 0:
            rowsO, pivO, MinvO, MO = gram_inverse('scalar', a, 0, sp.Rational(D1))
        else:
            rowsO, pivO, MinvO, MO = gram_inverse('cons', a, l)
        rowsF, pivF, MinvF, MF = gram_inverse('scalar', N - a, 0, s)
        assert len(pivF) == len(rowsF)
        iO = {r: i for i, r in enumerate(rowsO)}; iF = {r[0]: i for i, r in enumerate(rowsF)}
        for part in (sp.re, sp.im):
            Q = fmpq_mat(len(rowsO), len(rowsF))
            nz = False
            for (mx1, mz, mx2), c in entries.items():
                v = sp.Rational(part(c))
                if v != 0:
                    Q[iO[(mx1, mz)], iF[mx2]] = fmpq(int(v.p), int(v.q)); nz = True
            if not nz: continue
            # check that Q columns lie in range(M_O): solve on pivots and verify
            QP = fmpq_mat(len(pivO), len(rowsF))
            for i, p in enumerate(pivO):
                for j in range(len(rowsF)):
                    QP[i, j] = Q[p, j]
            U = MinvO*QP                       # coefficients on pivot columns
            # verify: M[:, piv] U == Q
            MOP = fmpq_mat(MO.nrows(), len(pivO))
            for i in range(MO.nrows()):
                for j, p in enumerate(pivO):
                    MOP[i, j] = MO[i, p]
            assert MOP*U == Q, "wavefunction not in the one-particle range (conservation violated?)"
            T = QP.transpose()*U*MinvF
            tr = sum((T[i, i] for i in range(T.nrows())), fmpq(0))
            tot += sp.Rational(int(tr.p), int(tr.q))
    return tot

# ---------------- structures at infinity ----------------
ZETA = (1, sp.I, 0)        # polarisation of X (null), zeta.zetabar = 2

def even_monomials(l, N, J):
    """T_n = (z1.x)^{l-n} (z1.zeta)^n (zeta.x)^{J-n} (x^2)^{p_n/2},  p_n = N - J - l + 2n (allowed iff p_n >= 0 even)"""
    zeta = sp.Matrix(ZETA)
    out = {}
    for n in range(l + 1):
        p = N - J - l + 2*n
        if p < 0 or p % 2 or J - n < 0: continue
        out[n] = sp.expand(dot(zv, xv)**(l-n)*dot(zv, zeta)**n*dot(zeta, xv)**(J-n)*dot(xv, xv)**(p//2))
    return out

def eps3(a, b, c):
    return sp.Matrix([a, b, c]).T.det() if False else sp.Matrix.hstack(a, b, c).det()

def odd_monomials(l, N, J):
    """eps(z1, zeta, x) * { (z1.x)^{l-1-n} (z1.zeta)^n (zeta.x)^{J-1-n} (x^2)^{p/2} },  n = 0..l-1"""
    zeta = sp.Matrix(ZETA)
    e = eps3(zv, zeta, xv)
    out = {}
    for n in range(l):
        p = N - 1 - (l - 1 - n) - (J - 1 - n)
        if p < 0 or p % 2 or J - 1 - n < 0: continue
        out[n] = sp.expand(e*dot(zv, xv)**(l-1-n)*dot(zv, zeta)**n*dot(zeta, xv)**(J-1-n)*dot(xv, xv)**(p//2))
    return out

def harmonic_z(q, l):
    if l < 2: return q
    lap = sum(sp.diff(q, Z[i], 2) for i in range(3))
    return sp.expand(q - dot(zv, zv)*lap/6)

def divergence(q, l):
    h = harmonic_z(q, l)
    return sp.expand(sum(sp.diff(h, X[i], Z[i]) for i in range(3)))

def conserved_combinations(mons, l):
    """null space of the divergence on span(mons); returns list of dicts n -> coefficient"""
    keys = sorted(mons)
    if not keys: return []
    cs = sp.symbols('c0:%d' % len(keys))
    div = sp.expand(sum(cs[i]*divergence(mons[k], l) for i, k in enumerate(keys)))
    if div == 0:
        return [{k: (1 if i == j else 0) for i, k in enumerate(keys)} for j in range(len(keys))]
    eqs = sp.Poly(div, *(X + Z)).coeffs()
    eqs = [sp.expand(e) for e in eqs]
    eqs = [e for e in eqs if e != 0]
    sol = sp.linsolve(eqs, cs)
    sol = list(sol)
    assert len(sol) == 1
    sol = sol[0]
    free = sorted(set().union(*[s.free_symbols for s in sol]) & set(cs), key=str)
    out = []
    for f in free:
        vec = [s.subs({g: (1 if g == f else 0) for g in free}) for s in sol]
        out.append({k: sp.nsimplify(vec[i]) for i, k in enumerate(keys)})
    return out
