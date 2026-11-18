"""
l4 (M6 step 2 piece): parity-even CONSERVED <T(x1,z1) T(x2,z2) O_{Delta,J}(x5,z5)> structures in d = 3 at general (complex Lorentzian)
position, as TPS-evaluable functions (for the crystal integrand), with exact Delta-dependent coefficients obtained numerically.

Basis (CPPR, embedding space, Lorentzian lightcone components (u, v, y), X = (1, x^2, x), Z = (0, 2 x.z, z)):
   V1^{a1} V2^{a2} V5^{a5} H12^{b12} H15^{b15} H25^{b25} / (X12^{(t1+t2-t5)/2} X15^{(t1+t5-t2)/2} X25^{(t2+t5-t1)/2}),  t_i = Delta_i + J_i,
   V1 = V_{1,25}, V2 = V_{2,51}, V5 = V_{5,12};  a1+b12+b15 = 2, a2+b12+b25 = 2, a5+b15+b25 = J.
Conservation of T1 (Delta1 = 3, degree 2 in z1):  C[f] = (3/2) d_x^mu d_{z mu} f - (1/2) z^nu d_{x nu} Box_z f = 0 (Todorov operator,
d = 3), evaluated with TPS in the 6 variables (x1, z1) at random points; same at point 2.  3d degeneracies (combinations vanishing
identically in d = 3) are quotiented out using the values of the structures at random points.
Checks: number of conserved, 1<->2 symmetric, parity-even structures (expected 2 for even J >= 2, generic Delta [Q: 3d <TTO> counting]),
and that the result is Delta-dependent (no Delta-independent solution).
Usage: python3 l4_tto_even.py -> cache/l4_tto_even.txt
"""
import os, sys, itertools, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l1_tps as T

DT = 3.0


def _is(a):
    return isinstance(a, np.ndarray) and a.shape == (T.N,)


def M(a, b):
    if _is(a) and _is(b):
        return T.mul(a, b)
    return a*b


def Dv(a, b):
    if _is(b):
        return T.div(a, b) if _is(a) else a*T.recip(b)
    return a/b


def P(a, p):
    if _is(a):
        return T.powr(a, p)
    return np.exp(p*np.log(a + 0j))


def ldot(a, b):
    return -0.5*(M(a[0], b[1]) + M(a[1], b[0])) + M(a[2], b[2])


_MODE = {'tps': False}


def _c(v):
    return T.const(v) if _MODE['tps'] else v


def EX(x):
    return (_c(1.0), ldot(x, x), x)


def EZ(x, z):
    return (_c(0.0), 2*ldot(x, z), z)


def tpsify(*vecs):
    """all components to TPS (numbers -> constants): scalar + TPS must never be added directly (numpy broadcast)."""
    return [[c if _is(c) else T.const(c) for c in v] for v in vecs]


def edot(A, B):
    return -0.5*(M(A[0], B[1]) + M(A[1], B[0])) + ldot(A[2], B[2])


def Xij(xi, xj):
    return -2*edot(EX(xi), EX(xj))


def Hs(xi, zi, xj, zj):
    Xi, Zi, Xj, Zj = EX(xi), EZ(xi, zi), EX(xj), EZ(xj, zj)
    return -2*(M(edot(Zi, Zj), edot(Xi, Xj)) - M(edot(Zi, Xj), edot(Zj, Xi)))


def Vs(xi, zi, xj, xk):
    Xi, Zi, Xj, Xk = EX(xi), EZ(xi, zi), EX(xj), EX(xk)
    return Dv(M(edot(Zi, Xj), edot(Xi, Xk)) - M(edot(Zi, Xk), edot(Xi, Xj)), edot(Xj, Xk))


def monomials(J):
    out = []
    for b12, b15, b25 in itertools.product(range(3), range(3), range(J + 1)):
        a1, a2, a5 = 2 - b12 - b15, 2 - b12 - b25, J - b15 - b25
        if min(a1, a2, a5) >= 0:
            out.append((a1, a2, a5, b12, b15, b25))
    return out


def structures(x1, z1, x2, z2, x5, z5, Delta, J, mons=None):
    """list of basis structure values (numbers or TPS)."""
    mons = mons or monomials(J)
    _MODE['tps'] = any(_is(c) for v in (x1, z1, x2, z2, x5, z5) for c in v)
    if _MODE['tps']:
        x1, z1, x2, z2, x5, z5 = tpsify(x1, z1, x2, z2, x5, z5)
    t1 = t2 = DT + 2; t5 = Delta + J
    V1, V2, V5 = Vs(x1, z1, x2, x5), Vs(x2, z2, x5, x1), Vs(x5, z5, x1, x2)
    H12, H15, H25 = Hs(x1, z1, x2, z2), Hs(x1, z1, x5, z5), Hs(x2, z2, x5, z5)
    den = M(M(P(Xij(x1, x2), -(t1 + t2 - t5)/2), P(Xij(x1, x5), -(t1 + t5 - t2)/2)), P(Xij(x2, x5), -(t2 + t5 - t1)/2))
    out = []
    for (a1, a2, a5, b12, b15, b25) in mons:
        v = den
        for base, k in ((V1, a1), (V2, a2), (V5, a5), (H12, b12), (H15, b15), (H25, b25)):
            for _ in range(k):
                v = M(v, base)
        out.append(v)
    return out


def _coef(a, e):
    return a[T._lookup(np.array([(np.array(e)*(T._B**np.arange(T.NV))).sum()]))[0]]


def conservation_row(pt, Delta, J, which=1):
    """C[f_k] at pt for every basis structure (conservation at point `which`), as a vector over (k, free z-component) rows."""
    x1, z1, x2, z2, x5, z5 = pt
    xs, zs_ = (x1, z1) if which == 1 else (x2, z2)
    X = [T.var(i, xs[i]) for i in range(3)]; Z = [T.var(3 + i, zs_[i]) for i in range(3)]
    args = (X, Z, x2, z2, x5, z5) if which == 1 else (x1, z1, X, Z, x5, z5)
    fs = structures(*args, Delta, J)
    ginv = np.array([[0, -2, 0], [-2, 0, 0], [0, 0, 1.0]])         # inverse lightcone metric
    def d(i):                                                       # exponent tuple for variable i
        e = [0]*T.NV; e[i] += 1; return e
    rows = []
    for f in fs:
        # first term: (3/2) g^{mn} d_{x m} d_{z n} f  (at base point: coefficient of x_m z_n)
        t1 = 0.0
        for m in range(3):
            for n in range(3):
                if ginv[m, n] != 0:
                    e = d(m); e[3 + n] += 1
                    t1 += ginv[m, n]*_coef(f, e)
        # second term: -(1/2) z^nu d_{x nu} Box_z f,  Box_z = g^{ab} d_{za} d_{zb}
        t2 = 0.0
        for nu_ in range(3):
            for a in range(3):
                for b in range(3):
                    if ginv[a, b] != 0:
                        e = d(nu_); e[3 + a] += 1; e[3 + b] += 1
                        mult = np.prod([math.factorial(k) for k in e])
                        t2 += zs_[nu_]*ginv[a, b]*_coef(f, e)*mult
        rows.append(1.5*t1 - 0.5*t2)
    return np.array(rows)


def rand_pt(rng):
    def null(rng):
        y = rng.normal(); u = rng.normal(); v = y*y/u                # -u v + y^2 = 0
        return np.array([u, v, y], complex)
    x = [rng.normal(size=3).astype(complex) for _ in range(3)]
    return (x[0], null(rng), x[1], null(rng), x[2], null(rng))


def conserved_even(Delta, J, npts=None, seed=3):
    rng = np.random.default_rng(seed)
    mons = monomials(J); K = len(mons)
    npts = npts or 3*K
    pts = [rand_pt(rng) for _ in range(npts)]
    Fv = np.array([structures(*p, Delta, J) for p in pts])                       # npts x K
    Cv = np.array([conservation_row(p, Delta, J, 1) for p in pts] + [conservation_row(p, Delta, J, 2) for p in pts])
    # 1<->2 symmetry: f(1<->2) = s f, s = (-1)^J for the standard basis [CPPR]; impose  f_sym - f = 0
    Sv = np.array([structures(p[2], p[3], p[0], p[1], p[4], p[5], Delta, J) for p in pts]) - Fv
    rn = lambda A_: A_/np.maximum(np.abs(A_).max(axis=1, keepdims=True), 1e-300)     # row scaling (null space unchanged)
    sc = np.maximum(np.abs(Fv).max(axis=1, keepdims=True), 1e-300)
    Fv, Cv, Sv = Fv/sc, rn(Cv), Sv/sc          # symmetry rows scaled like the structures (they may vanish identically)
    U, sF, Vh = np.linalg.svd(Fv); rF = int(np.sum(sF > 1e-9*sF[0]))
    B = Vh[:rF].conj().T                                   # K x rF: coordinates modulo 3d-vanishing combinations
    A = np.vstack([Cv @ B, Sv @ B])
    sA = np.linalg.svd(A, compute_uv=False)
    n0 = int(np.sum(sA < 1e-8*sA[0]))
    _, _, VhA = np.linalg.svd(A)
    sol = B @ VhA[len(VhA) - n0:].conj().T if n0 else np.zeros((K, 0))
    return dict(K=K, rank3d=rF, n_conserved_sym=n0, sv=sA[-4:]/sA[0], sol=sol, mons=mons)


if __name__ == '__main__':
    out = []
    def pr(s=""):
        print(s, flush=True); out.append(s)
    pr("# l4: parity-even conserved, 1<->2 symmetric <T T O_{Delta,J}> structures in d = 3 (numerical null space)")
    for J in (0, 2, 4, 6):
        for Dl in (7.3, 23.7, 101.1):
            r = conserved_even(Dl, J)
            pr(f"  J={J} Delta={Dl}: basis {r['K']}, independent in 3d {r['rank3d']}, conserved & symmetric: {r['n_conserved_sym']} "
               f"(smallest normalised singular values {np.array2string(r['sv'], precision=1)})")
    open(os.path.join(HERE, 'cache', 'l4_tto_even.txt'), 'w').write("\n".join(out) + "\n")
