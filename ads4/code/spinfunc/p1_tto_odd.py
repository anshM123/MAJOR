"""
p1 (parity-odd sector of route A', PARITY_ODD_ROUTEA.md sect. 1): parity-odd CONSERVED <T(x1,z1) T(x2,z2) O_{Delta,J}(x5,z5)> structures
in d = 3, covariant (embedding space, l4/l5 conventions), with exact Delta-dependent coefficients obtained numerically.

Conventions (as l4_tto_even): Lorentzian lightcone components (u, v, y), x.y = -(u v' + v u')/2 + y y'; embedding X = (X+, X-, x) =
(1, x^2, x), Z = (0, 2 x.z, z); A.B = -(A+ B- + A- B+)/2 + a.b;  X_ij = -2 X_i.X_j = (x_i - x_j)^2.  V_{i,jk}, H_ij as CLS-D / l4.
eps(A,B,C,D,E) := det of the 5 x 5 matrix of components (A+, A-, u, v, y) -- a REAL-coefficient polynomial, so a structure built with it
is real at real (Lorentzian, spacelike) configurations: its OPE coefficient is real in a unitary theory (Lorentzian-real basis).

Odd basis [CPPR 3d: every parity-odd structure is eps_ij = eps(Z_i, Z_j, X1, X2, X5) times a parity-even one]:
   type 12:  eps(Z1,Z2,X1,X2,X5) X12   x even monomial with spins (1, 1, J)
   type 15:  eps(Z1,Z5,X1,X2,X5) X15   x even monomial with spins (1, 2, J-1)
   type 25:  eps(Z2,Z5,X1,X2,X5) X25   x even monomial with spins (2, 1, J-1)
   all over the even denominator  X12^{(t1+t2-t5)/2} X15^{(t1+t5-t2)/2} X25^{(t2+t5-t1)/2}  times (X12 X15 X25)^{1/2}
   (eps_ij X_ij/sqrt(X12 X15 X25) has the homogeneity of an even (1,1,0)-type factor).
Conservation (Todorov operator, as l4), 1<->2 symmetry f(2<->1) = f, 3d relations quotiented numerically.
Expected: one conserved symmetric odd structure for J = 0, even J >= 2 and odd J >= 5; none for J = 1, 3 [flat: +- needs J >= 4].

Usage: python3 p1_tto_odd.py  -> cache/p1_tto_odd.txt  (counts at several Delta, J)
"""
import os, sys, itertools, math, functools, types
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np


# ------------------------------------------------------------------ embedding-level building blocks, generic over a namespace L
# L provides M (product), Dv (division), P (power), _is (is TPS), T (TPS module)
def edot(L, A, B):
    return -0.5*(L.M(A[0], B[1]) + L.M(A[1], B[0])) + (-0.5*(L.M(A[2], B[3]) + L.M(A[3], B[2])) + L.M(A[4], B[4]))


def emb_X(L, x):
    """x: (u, v, y) -> 5-vector (1, x^2, u, v, y)"""
    c1 = L.T.const(1.0) if any(L._is(c) for c in x) else 1.0
    x2 = -L.M(x[0], x[1]) + L.M(x[2], x[2])
    return [c1, x2, x[0], x[1], x[2]]


def emb_Z(L, x, z):
    tp = any(L._is(c) for c in list(x) + list(z))
    c0 = L.T.const(0.0) if tp else 0.0
    xz = -0.5*(L.M(x[0], z[1]) + L.M(x[1], z[0])) + L.M(x[2], z[2])
    zz = [c if (L._is(c) or not tp) else L.T.const(c) for c in z]
    return [c0, 2*xz, zz[0], zz[1], zz[2]]


def tpsify5(L, vs):
    if not any(L._is(c) for v in vs for c in v):
        return vs
    return [[c if L._is(c) else L.T.const(c) for c in v] for v in vs]


def eXX(L, Xi, Xj):
    return -2*edot(L, Xi, Xj)


def eV(L, Xi, Zi, Xj, Xk):
    num = L.M(edot(L, Zi, Xj), edot(L, Xi, Xk)) - L.M(edot(L, Zi, Xk), edot(L, Xi, Xj))
    return L.Dv(num, edot(L, Xj, Xk))


def eH(L, Xi, Zi, Xj, Zj):
    return -2*(L.M(edot(L, Zi, Zj), edot(L, Xi, Xj)) - L.M(edot(L, Zi, Xj), edot(L, Zj, Xi)))


_PERMS = None


def det5(L, rows):
    """5 x 5 determinant (Leibniz) of rows (lists of 5 numbers or TPS)."""
    global _PERMS
    if _PERMS is None:
        _PERMS = []
        for p in itertools.permutations(range(5)):
            inv = sum(1 for i in range(5) for j in range(i + 1, 5) if p[i] > p[j])
            _PERMS.append((p, -1 if inv % 2 else 1))
    tot = 0
    # expand along first row with 4x4 minors (cheaper): det = sum_j (-1)^j r0_j minor_0j
    for j in range(5):
        if not L._is(rows[0][j]) and rows[0][j] == 0:
            continue
        cols = [c for c in range(5) if c != j]
        sub = [[rows[i][c] for c in cols] for i in range(1, 5)]
        m4 = 0
        for p in itertools.permutations(range(4)):
            inv = sum(1 for a in range(4) for b in range(a + 1, 4) if p[a] > p[b])
            t = L.M(L.M(sub[0][p[0]], sub[1][p[1]]), L.M(sub[2][p[2]], sub[3][p[3]]))
            m4 = m4 + (-t if inv % 2 else t)
        tot = tot + ((-1)**j)*L.M(rows[0][j], m4)
    return tot


def odd_monomials(J):
    """(type, a1, a2, a5, b12, b15, b25): type in (12, 15, 25); spins used by eps: 12 -> (1,1,0), 15 -> (1,0,1), 25 -> (0,1,1)."""
    out = []
    for typ, (s1, s2, s5) in ((12, (1, 1, 0)), (15, (1, 0, 1)), (25, (0, 1, 1))):
        r1, r2, r5 = 2 - s1, 2 - s2, J - s5
        if r5 < 0:
            continue
        for b12, b15, b25 in itertools.product(range(3), range(3), range(J + 1)):
            a1, a2, a5 = r1 - b12 - b15, r2 - b12 - b25, r5 - b15 - b25
            if min(a1, a2, a5) >= 0:
                out.append((typ, a1, a2, a5, b12, b15, b25))
    return out


def structures_odd_emb(L, X1, Z1, X2, Z2, X5, Z5, Delta, J, mons=None, numerators=False):
    """values of the odd basis structures at embedding vectors (5-lists; numbers or TPS).  numerators=True: without the common
    denominator X12^{..} X15^{..} X25^{..} (X12 X15 X25)^{1/2}."""
    mons = mons or odd_monomials(J)
    X1, Z1, X2, Z2, X5, Z5 = tpsify5(L, [X1, Z1, X2, Z2, X5, Z5])
    DT = 3.0
    t1 = t2 = DT + 2; t5 = Delta + J
    V1, V2, V5 = eV(L, X1, Z1, X2, X5), eV(L, X2, Z2, X5, X1), eV(L, X5, Z5, X1, X2)
    H12, H15, H25 = eH(L, X1, Z1, X2, Z2), eH(L, X1, Z1, X5, Z5), eH(L, X2, Z2, X5, Z5)
    X12, X15, X25 = eXX(L, X1, X2), eXX(L, X1, X5), eXX(L, X2, X5)
    E = {12: L.M(det5(L, [Z1, Z2, X1, X2, X5]), X12), 15: L.M(det5(L, [Z1, Z5, X1, X2, X5]), X15),
         25: L.M(det5(L, [Z2, Z5, X1, X2, X5]), X25)}
    if numerators:
        den = L.T.const(1.0) if L._is(X12) else 1.0
    else:
        den = L.M(L.M(L.P(X12, -(t1 + t2 - t5)/2 - 0.5), L.P(X15, -(t1 + t5 - t2)/2 - 0.5)), L.P(X25, -(t2 + t5 - t1)/2 - 0.5))
    out = []
    pw = {}
    def p(nm, base, k):
        if (nm, k) not in pw:
            pw[(nm, k)] = (L.T.const(1.0) if L._is(base) else 1.0) if k == 0 else L.M(p(nm, base, k - 1), base)
        return pw[(nm, k)]
    for (typ, a1, a2, a5, b12, b15, b25) in mons:
        v = L.M(den, E[typ])
        for nm, base, k in (('V1', V1, a1), ('V2', V2, a2), ('V5', V5, a5), ('H12', H12, b12), ('H15', H15, b15), ('H25', H25, b25)):
            v = L.M(v, p(nm, base, k))
        out.append(v)
    return out


def structures_odd(L, x1, z1, x2, z2, x5, z5, Delta, J, mons=None, numerators=False):
    """physical points (lightcone 3-vectors)."""
    return structures_odd_emb(L, emb_X(L, x1), emb_Z(L, x1, z1), emb_X(L, x2), emb_Z(L, x2, z2), emb_X(L, x5), emb_Z(L, x5, z5),
                              Delta, J, mons, numerators)


# ------------------------------------------------------------------ conservation (Todorov operator, l4 code with the odd basis)
def conservation_row_odd(L, pt, Delta, J, which=1, mons=None):
    T = L.T
    x1, z1, x2, z2, x5, z5 = pt
    xs, zs_ = (x1, z1) if which == 1 else (x2, z2)
    X = [T.var(i, xs[i]) for i in range(3)]; Z = [T.var(3 + i, zs_[i]) for i in range(3)]
    args = (X, Z, x2, z2, x5, z5) if which == 1 else (x1, z1, X, Z, x5, z5)
    fs = structures_odd(L, *args, Delta, J, mons)
    ginv = np.array([[0, -2, 0], [-2, 0, 0], [0, 0, 1.0]])
    def d(i):
        e = [0]*T.NV; e[i] += 1; return e
    rows = []
    for f in fs:
        t1 = 0.0
        for m in range(3):
            for n in range(3):
                if ginv[m, n] != 0:
                    e = d(m); e[3 + n] += 1
                    t1 += ginv[m, n]*L._coef(f, e)
        t2 = 0.0
        for nu_ in range(3):
            for a in range(3):
                for b in range(3):
                    if ginv[a, b] != 0:
                        e = d(nu_); e[3 + a] += 1; e[3 + b] += 1
                        mult = np.prod([math.factorial(k) for k in e])
                        t2 += zs_[nu_]*ginv[a, b]*L._coef(f, e)*mult
        rows.append(1.5*t1 - 0.5*t2)
    return np.array(rows)


def conserved_odd_raw(L, Delta, J, npts=None, seed=3, sym=+1):
    rng = np.random.default_rng(seed)
    mons = odd_monomials(J); K = len(mons)
    if K == 0:
        return dict(K=0, rank3d=0, n_conserved_sym=0, sv=np.array([]), sol=np.zeros((0, 0)), mons=mons)
    npts = npts or 3*K
    pts = [L.rand_pt(rng) for _ in range(npts)]
    Fv = np.array([structures_odd(L, *p, Delta, J, mons) for p in pts])
    Cv = np.array([conservation_row_odd(L, p, Delta, J, 1, mons) for p in pts] + [conservation_row_odd(L, p, Delta, J, 2, mons) for p in pts])
    Sv = np.array([structures_odd(L, p[2], p[3], p[0], p[1], p[4], p[5], Delta, J, mons) for p in pts]) - sym*Fv
    rn = lambda A_: A_/np.maximum(np.abs(A_).max(axis=1, keepdims=True), 1e-300)
    sc = np.maximum(np.abs(Fv).max(axis=1, keepdims=True), 1e-300)
    Fv, Cv, Sv = Fv/sc, rn(Cv), Sv/sc
    U, sF, Vh = np.linalg.svd(Fv); rF = int(np.sum(sF > 1e-9*sF[0]))
    B = Vh[:rF].conj().T
    A = np.vstack([Cv @ B, Sv @ B])
    sA = np.linalg.svd(A, compute_uv=False)
    n0 = int(np.sum(sA < 1e-8*sA[0]))
    _, _, VhA = np.linalg.svd(A)
    sol = B @ VhA[len(VhA) - n0:].conj().T if n0 else np.zeros((K, 0))
    return dict(K=K, rank3d=rF, n_conserved_sym=n0, sv=sA[-4:]/sA[0], sol=sol, mons=mons, sF=sF[rF - 1:rF + 2]/sF[0])


@functools.lru_cache(maxsize=None)
def _Lred():
    """l4 helpers on the reduced TPS (6 vars, degree 3) with point-wise normalised powers (l5_spin._l4r)."""
    import l5_spin as S
    return S._l4r()


def expected_count(J):
    return 0 if J in (1, 3) else 1


@functools.lru_cache(maxsize=None)
def conserved_odd(Delta, J):
    """(mons, sol): sol K x 1, a REAL coefficient vector (normalised to unit max, sign: largest entry > 0) of the conserved,
    1<->2 symmetric parity-odd structure (None if there is none)."""
    L = _Lred()
    r = conserved_odd_raw(L, float(Delta), int(J))
    assert r['n_conserved_sym'] == expected_count(J), (Delta, J, r['n_conserved_sym'], r['sv'])
    if r['n_conserved_sym'] == 0:
        return tuple(r['mons']), None, r['sv']
    sol = r['sol'][:, 0]
    ph = sol[np.argmax(np.abs(sol))]
    sol = sol/ph
    assert np.abs(sol.imag).max() < 1e-7, ('not real', Delta, J, np.abs(sol.imag).max())
    return tuple(r['mons']), sol.real[:, None], r['sv']


if __name__ == '__main__':
    out = []
    def pr(s=""):
        print(s, flush=True); out.append(s)
    L = _Lred()
    pr("# p1: parity-odd conserved <T T O_{Delta,J}> structures in d = 3 (numerical null space); expected 1 (J = 0, 2, 4, 5, 6, ...), 0 (J = 1, 3)")
    for J in (0, 1, 2, 3, 4, 5, 6, 7, 9):
        for Dl in (7.3, 23.7, 101.1, 1003.7):
            r = conserved_odd_raw(L, Dl, J)
            ra = conserved_odd_raw(L, Dl, J, sym=-1)
            pr(f"  J={J} Delta={Dl}: odd basis {r['K']}, independent in 3d {r['rank3d']}, conserved & symmetric {r['n_conserved_sym']} "
               f"(antisymmetric {ra['n_conserved_sym']}); smallest sv {np.array2string(r['sv'], precision=1)}")
    open(os.path.join(HERE, 'cache', 'p1_tto_odd.txt'), 'w').write("\n".join(out) + "\n")
