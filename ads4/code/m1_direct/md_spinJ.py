"""md_spinJ.py -- point-wise EE k = 2 CLS-D integrand on D' for a block of spin J (even), ALL conserved parity-even <T T O_{Delta,J}>
structures, own numpy implementation of the spin contractions (null-circle Fourier encoding, unit-tested against explicit tensors).

  Phi_EE2[a, c] = Phi_scal(J=0) * u'v'/(m12^2 m34^2) * |m12 m34/(m15 m25 m35 m45)|^{J/2} * prefh_J[a, c],
  prefh_J[a, c] = < F_a(z2, z5) , Ehat(z2, z4)^2 , G_c(z4, z5) >   (spin-2 contractions at x2, x4; spin-J contraction at x5 = e),
  F_a = numerator of structure a of <T(x1, z1 = x10) T(x2, z2) O_{Delta,J}(e, z5)>, G_c = same for <T(x3, z3 = x30) T(x4, z4) O_{3-Delta,J}(e, z5)>,
  Ehat = (H24 - 2 V_{2,04} V_{4,02})/(X20 X40)  (the EE dual structure divided by V042^2).
The J-dependent m-powers are the 3-point denominators X12^{-(10-Delta-J)/2} X15^{-(Delta+J)/2} ... minus their J = 0 values
(w0_ee_common._gs: -((10-J)/2) log m12 - (J/2) log(m15 m25) - ((7-J)/2) log m34 - ((3+J)/2) log(m35 m45)).
Structures: monomial basis and conserved null space from spinfunc/l4_tto_even (imported read-only; conserved_even writes nothing).
V, H: CPPR embedding formulas (identical to repro_engine Vf/Hf, checked).
At J = 0 this reproduces rd_core's Q/R (md_dquad.Q2) up to a constant (unit test 'j0').
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
import numpy as np
from math import comb

GL = np.array([[0.0, -0.5, 0.0], [-0.5, 0.0, 0.0], [0.0, 0.0, 1.0]])   # lightcone metric (u, v, y)
GI = np.linalg.inv(GL)


def dot(p, q):
    return -0.5*(p[0]*q[1] + p[1]*q[0]) + p[2]*q[2]


def sub(p, q):
    return [p[i] - q[i] for i in range(3)]


def sq(p, q):
    r = sub(p, q); return dot(r, r)


def Vf(xi, z, xj, xk):
    return (dot(z, sub(xj, xi))*sq(xi, xk) - dot(z, sub(xk, xi))*sq(xi, xj))/sq(xj, xk)


def Hf(xi, zi, xj, zj):
    return dot(zi, zj)*sq(xi, xj) + 2*dot(zi, sub(xj, xi))*dot(zj, sub(xi, xj))


def null_circle(l):
    """2l+1 equispaced null vectors w(phi) = (t, x, y) = (1, cos, sin) in lightcone comps (u, v, y) = (1+cos, 1-cos, sin)."""
    ph = 2*np.pi*np.arange(2*l + 1)/(2*l + 1)
    return ph, np.stack([1 + np.cos(ph), 1 - np.cos(ph), np.sin(ph)], axis=1)       # (2l+1, 3)


def kappa(l):
    n = np.arange(-l, l + 1)
    return np.array([(-1.0)**(l + k)*2.0**(-l)*comb(2*l, l + k) for k in n])


def fourier(vals, axis, l):
    """samples on null_circle(l) along `axis` -> Fourier coefficients c_n, n = -l..l (same axis)."""
    ph, _ = null_circle(l)
    n = np.arange(-l, l + 1)
    Fm = np.exp(-1j*np.outer(n, ph))/(2*l + 1)                      # (2l+1 modes, 2l+1 samples)
    v = np.moveaxis(vals, axis, 0)
    out = np.tensordot(Fm, v, axes=([1], [0]))
    return np.moveaxis(out, 0, axis)


def contract(fa, fb, l, axa, axb):
    """<A, B> over one spin-l index: A, B given by Fourier coefficients along axa, axb: sum_n A_n B_{-n}/kappa_n."""
    k = kappa(l)
    A = np.moveaxis(fa, axa, -1); B = np.moveaxis(fb, axb, -1)[..., ::-1]
    return A, B/k


# ------------------------------------------------------------------ structures
_CACHE = {}


def conserved(Delta, J):
    key = (round(float(Delta), 12), int(J))
    if key not in _CACHE:
        import l4_tto_even as L4
        r = L4.conserved_even(float(Delta), int(J))
        sol = r['sol']
        U, sv, _ = np.linalg.svd(np.hstack([sol.real, sol.imag]), full_matrices=False)
        nr = int(np.sum(sv > 1e-8*sv[0]))
        sol = U[:, :nr]
        sol = sol/np.abs(sol).max(axis=0, keepdims=True)
        _CACHE[key] = (list(r['mons']), sol, r['n_conserved_sym'])
    return _CACHE[key]


def tto_num(xa, za, xb, zbs, x5, z5s, J, mons, sol):
    """numerators of the conserved structures (columns of sol) at points (arrays N): za fixed vector field, zb in zbs (nb, 3),
    z5 in z5s (n5, 3) -> (nstruct, nb, n5, N)."""
    def lin(f, Z):
        b = [np.broadcast_to(f([1.0 if i == m else 0.0 for i in range(3)]), np.shape(xa[0])) for m in range(3)]
        return Z[:, 0, None]*b[0][None] + Z[:, 1, None]*b[1][None] + Z[:, 2, None]*b[2][None]
    V1 = Vf(xa, za, xb, x5)                                                       # (N)
    V2 = lin(lambda B: Vf(xb, B, x5, xa), zbs)                                    # (nb, N)
    V5 = lin(lambda B: Vf(x5, B, xa, xb), z5s)                                    # (n5, N)
    H12 = lin(lambda B: Hf(xa, za, xb, B), zbs)                                   # (nb, N)
    H15 = lin(lambda B: Hf(xa, za, x5, B), z5s)                                   # (n5, N)
    Hb = np.stack([np.stack([Hf(xb, [1.0*(i == m) for i in range(3)], x5, [1.0*(i == k) for i in range(3)]) for k in range(3)])
                   for m in range(3)])                                            # (3, 3, N)
    H25 = np.einsum('im,jk,mkn->ijn', zbs, z5s, Hb)                               # (nb, n5, N)
    nb, n5, N = zbs.shape[0], z5s.shape[0], np.shape(V1)[-1]
    res = np.zeros((sol.shape[1], nb, n5, N))
    pw = lambda b, k: [np.ones_like(b)] + [b**j for j in range(1, k + 1)]
    V1p, V2p, H12p, V5p, H15p, H25p = pw(V1, 2), pw(V2, 2), pw(H12, 2), pw(V5, J), pw(H15, 2), pw(H25, 2)
    for i, (a1, a2, a5, b12, b15, b25) in enumerate(mons):
        if np.all(np.abs(sol[i]) < 1e-14):
            continue
        A = V1p[a1][None]*V2p[a2]*H12p[b12]                       # (nb, N)
        Bv = V5p[a5]*H15p[b15]                                    # (n5, N)
        Mv = A[:, None]*Bv[None, :]*H25p[b25]
        for s_ in range(sol.shape[1]):
            if sol[i, s_] != 0:
                res[s_] += sol[i, s_]*Mv
    return res


def prefh(g, Delta, J):
    """prefh_J[a, c] (nsa, nsc, N) at geometry values g (rd_core.geometry dict, arrays)."""
    re = lambda v: [np.real(c) if np.ndim(c) else float(np.real(c)) for c in v]     # all positions are real on D'
    x1, x2, x3, x4, x0, e = (re(g[k]) for k in ('x1', 'x2', 'x3', 'x4', 'x0', 'e'))
    x10, x30 = sub(x1, x0), sub(x3, x0)
    ph2, W2 = null_circle(2); phJ, WJ = null_circle(J)
    monsD, solD, _ = conserved(Delta, J)
    monsT, solT, _ = conserved(3.0 - Delta, J)
    F = tto_num(x1, x10, x2, W2, e, WJ, J, monsD, solD)                           # (na, 5, 2J+1, N)
    G = tto_num(x3, x30, x4, W2, e, WJ, J, monsT, solT)                           # (nc, 5, 2J+1, N)
    # Ehat(z2, z4) on the null circles, squared
    inv = 1.0/(sq(x2, x0)*sq(x4, x0))
    B3 = [[1.0*(i == m) for i in range(3)] for m in range(3)]
    H24 = np.stack([np.stack([Hf(x2, B3[m], x4, B3[k]) for k in range(3)]) for m in range(3)])
    V204 = np.stack([Vf(x2, B3[m], x0, x4) for m in range(3)]); V402 = np.stack([Vf(x4, B3[m], x0, x2) for m in range(3)])
    Em = inv*(H24 - 2*np.einsum('mn,kn->mkn', V204, V402))                       # (3, 3, N)
    Es = np.einsum('im,jk,mkn->ijn', W2, W2, Em)                                  # (5, 5, N)
    D = Es*Es
    # Fourier transforms on every null-circle axis
    Fh = fourier(fourier(F, 1, 2), 2, J)
    Gh = fourier(fourier(G, 1, 2), 2, J)
    Dh = fourier(fourier(D, 0, 2), 1, 2)
    k2, kJ = kappa(2), kappa(J)
    # <F, D>_z2 : sum_n F_n D_{-n, .}/k_n
    FD = np.einsum('aipn,iqn->apqn', Fh, Dh[::-1]/k2[:, None, None])               # (na, 2J+1, 5[z4 modes], N)
    # <., G>_z4 and <.,.>_z5
    out = np.einsum('apqn,cqrn->acprn', FD, Gh[:, ::-1]/k2[None, :, None, None])  # (na, nc, 2J+1, 2J+1, N)
    out = np.einsum('acppn->acn', out[:, :, :, ::-1]/kJ[None, None, None, :, None]) if False else \
        np.einsum('acpn->acn', np.stack([out[:, :, p, 2*J - p, :]/kJ[2*J - p] for p in range(2*J + 1)], axis=2))
    return out


def jfactor(g, J):
    up = g['m12']*g['m34']/(g['m13']*g['m24']); vp = g['m14']*g['m23']/(g['m13']*g['m24'])
    f = up*vp/(g['m12']**2*g['m34']**2)
    if J:
        f = f*np.abs(g['m12']*g['m34']/(g['m15']*g['m25']*g['m35']*g['m45']))**(J/2)
    return f


def ratio(g, Delta, J):
    """Phi_EE2[a, c]/Phi_scal (nsa, nsc, N), real part (Phi is real on D')."""
    return np.real(prefh(g, Delta, J)*jfactor(g, J)[None, None])


# ------------------------------------------------------------------ unit tests
def _test_contract(rng):
    """<T, S> for S = (w'.w)^l, w' a random null vector (S = w'^{(x)l} is traceless): must equal T(w') for any T; T = prod (a_i . w)."""
    out = []
    for l in (2, 4, 8):
        A = [rng.normal(size=3) for _ in range(l)]
        th = rng.uniform(0, 2*np.pi); lam = rng.uniform(0.5, 2)
        wp = lam*np.array([1 + np.cos(th), 1 - np.cos(th), np.sin(th)])
        ph, W = null_circle(l)
        fa = np.prod([np.array([dot(w, a) for w in W]) for a in A], axis=0)
        fs = np.array([dot(w, wp) for w in W])**l
        val = np.sum(fourier(fa, 0, l)*fourier(fs, 0, l)[::-1]/kappa(l))
        ex = np.prod([dot(wp, a) for a in A])
        out.append((l, val, ex))
    return out


if __name__ == '__main__':
    rng = np.random.default_rng(5)
    for l, v, ex in _test_contract(rng):
        print('contraction test l =', l, v, ex)


# ------------------------------------------------------------------ k = 3 piece (C^_3 = Psi_{3,nu}[o_H q], j' = 1, f = -1)
def _emb_X(x):
    return np.stack([np.ones_like(x[0]), dot(x, x), x[0], x[1], x[2]], axis=-1)


def _emb_Z(x, z):
    sh = np.shape(x[0])
    z = [np.broadcast_to(np.asarray(c, float), sh) for c in z]
    return np.stack([np.zeros(sh), 2*dot(x, z), z[0], z[1], z[2]], axis=-1)


def _cof(R1, R2, R3, R4):
    """cofactor vectors c (N, 5): det([r0; R1; R2; R3; R4]) = r0 . c  (rows (N, 5))."""
    M = np.stack([np.zeros_like(R1), R1, R2, R3, R4], axis=1)          # (N, 5, 5)
    c = []
    for i in range(5):
        Mi = M.copy(); Mi[:, 0, :] = 0.0; Mi[:, 0, i] = 1.0
        c.append(np.linalg.det(Mi))
    return np.stack(c, axis=-1)


def prefh3(g, Delta, J):
    """< F_a, q (x) o_H, G_c > / v042^4  (nsa, nsc, N); l9_points.build_spin kind 'C3H' re-implemented point-wise."""
    re = lambda v: [np.real(c) if np.ndim(c) else float(np.real(c)) for c in v]
    x1, x2, x3, x4, x0, e, z0 = (re(g[k]) for k in ('x1', 'x2', 'x3', 'x4', 'x0', 'e', 'z0'))
    N = np.shape(x1[0])[-1]
    z0 = [np.broadcast_to(np.asarray(c, float), (N,)) for c in z0]
    x10, x30 = sub(x1, x0), sub(x3, x0)
    ph2, W2 = null_circle(2); phJ, WJ = null_circle(J)
    monsD, solD, _ = conserved(Delta, J)
    monsT, solT, _ = conserved(3.0 - Delta, J)
    F = tto_num(x1, x10, x2, W2, e, WJ, J, monsD, solD)
    G = tto_num(x3, x30, x4, W2, e, WJ, J, monsT, solT)
    B3 = [[1.0*(i == m) for i in range(3)] for m in range(3)]
    lin = lambda f: np.stack([f(B3[m]) for m in range(3)])                         # (3, N)
    H20 = W2 @ lin(lambda B: Hf(x2, B, x0, z0)); H40 = W2 @ lin(lambda B: Hf(x4, B, x0, z0))   # (5, N)
    H24b = np.stack([np.stack([Hf(x2, B3[m], x4, B3[k]) for k in range(3)]) for m in range(3)])
    H24 = np.einsum('im,jk,mkn->ijn', W2, W2, H24b)
    X20, X40, X24 = sq(x2, x0), sq(x4, x0), sq(x2, x4)
    v042 = Vf(x0, z0, x4, x2)
    q = (2*H20[:, None]*H40[None, :] + v042**2*H24)/(X20*X40)
    # o_H = eps(Z2,X2,X4,X0,Z0) H40/sqrt(X24 X20 X40^3) + eps(Z4,X4,X2,X0,Z0) H20/sqrt(X24 X40 X20^3)
    EX2, EX4, EX0, EZ0 = _emb_X(x2), _emb_X(x4), _emb_X(x0), _emb_Z(x0, z0)
    c = _cof(EX2, EX4, EX0, EZ0)                                                   # eps(Z2, X2, X4, X0, Z0) = Z2 . c
    e2 = np.stack([np.einsum('ni,ni->n', _emb_Z(x2, [W2[j, 0], W2[j, 1], W2[j, 2]]), c) for j in range(5)])
    e4 = np.stack([-np.einsum('ni,ni->n', _emb_Z(x4, [W2[j, 0], W2[j, 1], W2[j, 2]]), c) for j in range(5)])
    n2 = 1.0/np.sqrt(np.abs(X24*X20*X40**3)); n4 = 1.0/np.sqrt(np.abs(X24*X40*X20**3))
    o = (e2*n2)[:, None]*H40[None, :] + H20[:, None]*(e4*n4)[None, :]
    D = q*o/v042**4
    Fh = fourier(fourier(F, 1, 2), 2, J); Gh = fourier(fourier(G, 1, 2), 2, J); Dh = fourier(fourier(D, 0, 2), 1, 2)
    k2, kJ = kappa(2), kappa(J)
    FD = np.einsum('aipn,iqn->apqn', Fh, Dh[::-1]/k2[:, None, None])
    out = np.einsum('apqn,cqrn->acprn', FD, Gh[:, ::-1]/k2[None, :, None, None])
    return np.einsum('acpn->acn', np.stack([out[:, :, p, 2*J - p, :]/kJ[2*J - p] for p in range(2*J + 1)], axis=2))


def jsign(g):
    """the j' = 1 celestial sign: sign det[z0; x10; x30] (3-d, lightcone components; a global sign convention is irrelevant)."""
    re = lambda v: [np.real(c) if np.ndim(c) else float(np.real(c)) for c in v]
    x1, x3, x0, z0 = (re(g[k]) for k in ('x1', 'x3', 'x0', 'z0'))
    a, b = sub(x1, x0), sub(x3, x0)
    N = np.shape(a[0])[-1]
    z = [np.broadcast_to(np.asarray(c, float), (N,)) for c in z0]
    det = z[0]*(a[1]*b[2] - a[2]*b[1]) - z[1]*(a[0]*b[2] - a[2]*b[0]) + z[2]*(a[0]*b[1] - a[1]*b[0])
    return np.sign(det)


def ratio3(g, Delta, J):
    """Phi_C3H[a, c]/Phi_scal: [u'v'/(v'-u')] m12^-2 m34^-2 |m12 m34/(m15 m25 m35 m45)|^{J/2} prefh3 * sign_{j'=1}."""
    up = g['m12']*g['m34']/(g['m13']*g['m24']); vp = g['m14']*g['m23']/(g['m13']*g['m24'])
    f = up*vp/(vp - up)/(g['m12']**2*g['m34']**2)
    if J:
        f = f*np.abs(g['m12']*g['m34']/(g['m15']*g['m25']*g['m35']*g['m45']))**(J/2)
    return np.real(prefh3(g, Delta, J)*(f*jsign(g))[None, None])
