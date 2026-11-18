"""
l5 (M6 step 2, round 3): machinery for the spinning (<TTTT>) CLS-D crystal integrand, parity-even sector.

* bm: batched TPS product (l1_tps engine, arrays of shape (..., N)).
* Embedding-space V, H (CLS-D eq. CFT_building_blocks_VH, same conventions as l4_tto_even) on batched TPS, linear in the polarisations,
  so that values at many sampled polarisations come from 3 basis evaluations.
* Index contraction of symmetric traceless spin-l tensors in d = 3 (Lorentzian metric) from values on the null circle
  w(phi) = (t, x, y) = (1, cos phi, sin phi):  <T, S> = sum_{n=-l..l} T_n S_{-n} / kappa_n,  T_n = Fourier coefficients of T(w(phi)),
  kappa_n = (-1)^{l+n} 2^{-l} C(2l, l+n)  [D: T = sum_k c_k w_k^{(x)l}, (w.w')^l = (cos(phi-phi') - 1)^l = sum kappa_n e^{in(phi-phi')}].
* Conserved parity-even <T T O_{Delta,J}> coefficient vectors (l4_tto_even.conserved_even) computed on a reduced TPS (6 variables, degree 3:
  enough for the Todorov operator) with the denominators normalised point-wise (no overflow at |Delta| ~ 10^3; a per-point constant does not
  change the null space).
Self-tests: python3 l5_spin.py  (contraction identities; scalar Q_J reproduced by the Fourier contraction on the crystal TPS; batched V/H
vs l4 at a numeric point).
"""
import os, sys, types, math, functools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import scipy.sparse as sp
from scipy.special import comb
import l1_tps as T

N = T.N
_PM = sp.csr_matrix((np.ones(len(T._K)), (np.arange(len(T._K)), T._K)), shape=(len(T._K), N))
_PMT = _PM.T.tocsr()


def bm(a, b):
    """batched TPS product; a, b arrays (..., N) (broadcast)."""
    a = np.asarray(a, complex); b = np.asarray(b, complex)
    sh = np.broadcast_shapes(a.shape[:-1], b.shape[:-1])
    a2 = np.broadcast_to(a, sh + (N,)).reshape(-1, N); b2 = np.broadcast_to(b, sh + (N,)).reshape(-1, N)
    out = np.empty((a2.shape[0], N), complex)
    ch = max(1, int(1.5e7//len(T._I)))
    for s in range(0, a2.shape[0], ch):
        w = a2[s:s + ch, T._I]*b2[s:s + ch, T._J]
        out[s:s + ch] = (_PMT @ w.T).T
    return out.reshape(sh + (N,))


def C(c):
    return T.const(c)


def tv(v):
    """numeric 3-vector -> list of const TPS."""
    return [C(c) for c in v]


def ld(a, b):
    return -0.5*(bm(a[0], b[1]) + bm(a[1], b[0])) + bm(a[2], b[2])


def EX(x):
    return (C(1.0), ld(x, x), x)


def ed_XX(A, B):
    return -0.5*(B[1] + A[1]) + ld(A[2], B[2])            # both have first component 1


def ed_ZX(Z, X):                                          # Z = (0, 2x.z, z), X = (1, x^2, x)
    return -0.5*Z[1] + ld(Z[2], X[2])


def ed_ZZ(Z1, Z2):
    return ld(Z1[2], Z2[2])


def EZ(x, z):
    return (None, 2*ld(x, z), z)


def Vfun(xi, zi, xj, xk):
    """V_{i,jk} = (Z_i.X_j X_i.X_k - Z_i.X_k X_i.X_j)/(X_j.X_k); zi may be batched."""
    Xi, Xj, Xk, Zi = EX(xi), EX(xj), EX(xk), EZ(xi, zi)
    num = bm(ed_ZX(Zi, Xj), ed_XX(Xi, Xk)) - bm(ed_ZX(Zi, Xk), ed_XX(Xi, Xj))
    return bm(num, T.recip(ed_XX(Xj, Xk)))


def Hfun(xi, zi, xj, zj):
    """H_ij = -2(Z_i.Z_j X_i.X_j - Z_i.X_j Z_j.X_i)."""
    Xi, Xj, Zi, Zj = EX(xi), EX(xj), EZ(xi, zi), EZ(xj, zj)
    return -2*(bm(ed_ZZ(Zi, Zj), ed_XX(Xi, Xj)) - bm(ed_ZX(Zi, Xj), ed_ZX(Zj, Xi)))


def Xfun(xi, xj):
    return -2*ed_XX(EX(xi), EX(xj))


BASIS = [tv(np.eye(3)[m]) for m in range(3)]             # lightcone basis vectors (u, v, y)


def null_circle(l):
    """2l+1 null vectors w(phi_k) in lightcone components (u, v, y) = (1 - cos, 1 + cos, sin)."""
    n = 2*l + 1; ph = 2*np.pi*np.arange(n)/n
    return np.stack([1 - np.cos(ph), 1 + np.cos(ph), np.sin(ph)], axis=1)


def kappa(l):
    n = np.arange(-l, l + 1)
    return (-1.0)**(l + n)*2.0**(-l)*comb(2*l, l + n)


def fourier(vals, axis, l):
    """Fourier coefficients n = -l..l (in that order) of samples on null_circle(l) along `axis`."""
    nn = 2*l + 1
    F = np.fft.fft(vals, axis=axis)/nn
    return np.take(F, np.arange(-l, l + 1) % nn, axis=axis)


def lin_samples(fbasis, W):
    """fbasis: list of 3 arrays (values for z = basis vectors, linear in z); W: (n, 3) sample vectors -> (n, ...) array."""
    fb = np.stack(fbasis, axis=0)
    return np.tensordot(W, fb, axes=([1], [0]))


# ------------------------------------------------------------------ conserved <TTO> coefficients on a reduced TPS
@functools.lru_cache(maxsize=None)
def _l4r():
    src = open(os.path.join(HERE, 'l1_tps.py')).read().replace("NV, DEG = 7, 6", "NV, DEG = 6, 3")
    T3 = types.ModuleType('l1_tps_63'); exec(compile(src, 'l1_tps_63', 'exec'), T3.__dict__)
    saved = sys.modules.get('l1_tps'); sys.modules['l1_tps'] = T3
    try:
        L = types.ModuleType('l4r'); L.__dict__['__file__'] = os.path.join(HERE, 'l4_tto_even.py')
        exec(compile(open(os.path.join(HERE, 'l4_tto_even.py')).read(), 'l4r', 'exec'), L.__dict__)
    finally:
        sys.modules['l1_tps'] = saved
    def Pn(a, p):                                         # point-wise normalised power (constant factor per point)
        if L._is(a):
            return T3.powr(a/abs(a[0]), p)
        return np.exp(p*np.log(a/abs(a) + 0j))
    L.P = Pn
    return L


@functools.lru_cache(maxsize=None)
def conserved(Delta, J):
    """(mons, sol): sol is K x n (n = 1 for J = 0, 2; 2 for even J >= 4); columns normalised to unit max."""
    L = _l4r()
    r = L.conserved_even(float(Delta), int(J))
    sol = r['sol']
    # the null space is spanned by real vectors (real Delta): rotate the (complex, SVD-arbitrary) basis to a real orthonormal one
    U, sv_, _ = np.linalg.svd(np.hstack([sol.real, sol.imag]), full_matrices=False)
    nreal = int(np.sum(sv_ > 1e-8*sv_[0]))
    assert nreal == sol.shape[1], ('null space not real', Delta, J, sv_)
    sol = U[:, :nreal]
    assert r['n_conserved_sym'] == (1 if J < 4 else 2), (Delta, J, r['n_conserved_sym'], r['sv'])
    sol = sol/np.abs(sol).max(axis=0, keepdims=True)
    return tuple(r['mons']), sol, r['sv']


def tto_samples(xa, za, xb, zbs, x5, z5s, J, mons, sol):
    """parity-even <T(xa,za) T(xb,zb) O(x5,z5)> NUMERATORS (all monomials of l4 times sol; the common denominator
    X_ab^{-(10-Delta-J)/2} X_a5^{-(Delta+J)/2} X_b5^{-(Delta+J)/2} is NOT included) for za fixed, zb in zbs (nb x 3), z5 in z5s (n5 x 3).
    Returns array (nstruct, nb, n5, N)."""
    V1 = Vfun(xa, za, xb, x5)                                                   # (N)
    V2 = lin_samples([Vfun(xb, BASIS[m], x5, xa) for m in range(3)], zbs)       # (nb, N)
    V5 = lin_samples([Vfun(x5, BASIS[m], xa, xb) for m in range(3)], z5s)       # (n5, N)
    H12 = lin_samples([Hfun(xa, za, xb, BASIS[m]) for m in range(3)], zbs)      # (nb, N)
    H15 = lin_samples([Hfun(xa, za, x5, BASIS[m]) for m in range(3)], z5s)      # (n5, N)
    H25b = [[Hfun(xb, BASIS[m], x5, BASIS[k]) for k in range(3)] for m in range(3)]
    H25 = np.einsum('im,jk,mkn->ijn', zbs, z5s, np.array(H25b))                 # (nb, n5, N)
    def pw(base, kmax):
        out = [np.broadcast_to(C(1.0), base.shape).copy()]
        for _ in range(kmax):
            out.append(bm(out[-1], base))
        return out
    V1p, V2p, V5p, H12p, H15p, H25p = pw(V1, 2), pw(V2, 2), pw(V5, J), pw(H12, 2), pw(H15, 2), pw(H25, 2)
    nb, n5 = zbs.shape[0], z5s.shape[0]
    res = np.zeros((sol.shape[1], nb, n5, N), complex)
    for i, (a1, a2, a5, b12, b15, b25) in enumerate(mons):
        if np.all(np.abs(sol[i]) < 1e-14):
            continue
        A = bm(bm(V1p[a1], V2p[a2]), H12p[b12])                                # (nb, N)
        B = bm(V5p[a5], H15p[b15])                                              # (n5, N)
        Mv = bm(bm(A[:, None, :], B[None, :, :]), H25p[b25])                    # (nb, n5, N)
        res += sol[i][:, None, None, None]*Mv[None]
    return res


# ------------------------------------------------------------------ self-tests
def _tests():
    out = []
    def pr(s=""):
        print(s, flush=True); out.append(s)
    rng = np.random.default_rng(1)
    def rnull():
        y = rng.normal() + 1j*rng.normal(); u = rng.normal() + 1j*rng.normal(); return np.array([u, y*y/u, y])
    dl = lambda a, b: -0.5*(a[0]*b[1] + a[1]*b[0]) + a[2]*b[2]
    pr("# l5 self-tests")
    for l in (1, 2, 3, 6, 12, 20, 28):
        a, b = rnull(), rnull()
        W = null_circle(l)
        fa = np.array([dl(w, a)**l for w in W]); fb = np.array([dl(w, b)**l for w in W])
        Fa, Fb = fourier(fa, 0, l), fourier(fb, 0, l)
        val = np.sum(Fa*Fb[::-1]/kappa(l))
        pr(f"  contraction <(a.z)^l,(b.z)^l> vs (a.b)^l, l={l}: rel err {abs(val/dl(a, b)**l - 1):.1e}")
    # l = 2 with a non-null-built tensor: T(z) = (a.z)(b.z) (traceless part) vs S = (c.z)^2:  <T,S> = (a.c)(b.c) (c null)
    a, b, c = rng.normal(size=3) + 0j, rng.normal(size=3) + 0j, rnull()
    W = null_circle(2)
    val = np.sum(fourier(np.array([dl(w, a)*dl(w, b) for w in W]), 0, 2)*fourier(np.array([dl(w, c)**2 for w in W]), 0, 2)[::-1]/kappa(2))
    pr(f"  l=2 non-null T: rel err {abs(val/(dl(a, c)*dl(b, c)) - 1):.1e}")
    # batched V/H vs l4 at a numeric point
    L = _l4r()
    x1, x2, x5 = [rng.normal(size=3) + 0j for _ in range(3)]; z1, z2, z5 = rnull(), rnull(), rnull()
    import l4_tto_even as L4
    for nm, mine, ref in (('V125', Vfun(tv(x1), tv(z1), tv(x2), tv(x5))[0], L4.Vs(x1, z1, x2, x5)),
                          ('H25', Hfun(tv(x2), tv(z2), tv(x5), tv(z5))[0], L4.Hs(x2, z2, x5, z5))):
        pr(f"  {nm}: batched TPS vs l4: rel err {abs(mine/ref - 1):.1e}")
    mons, sol, sv = conserved(101.1, 4)
    st = tto_samples(tv(x1), tv(z1), tv(x2), z2[None, :], tv(x5), z5[None, :], 4, mons, sol)[:, 0, 0, 0]
    L4.P = lambda a, p: 1.0                                                   # numerators only
    ref = np.array(L4.structures(x1, z1, x2, z2, x5, z5, 101.1, 4)) @ sol
    pr(f"  <TTO> numerators (J=4, two conserved structures) vs l4: rel err {np.max(np.abs(st/ref - 1)):.1e}; null-space sv {sv}")
    return out


if __name__ == '__main__':
    o = _tests()
    os.makedirs(os.path.join(HERE, 'cache'), exist_ok=True)
    open(os.path.join(HERE, 'cache', 'l5_spin_tests.txt'), 'w').write("\n".join(o) + "\n")
