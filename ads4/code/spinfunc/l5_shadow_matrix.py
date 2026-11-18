"""
l5_shadow_matrix (M6 step 2, item 1 of the 20.9 recipe): the spinning shadow-coefficient MATRIX beta_{cb}(Delta, J) of <T T O_{Delta,J}>,
parity-even sector, from the Laplace expansion of the shadow integral around x5 = e to relative O(1/Delta^2) (l3 route, validated there on
the scalar beta_{Delta,J}).

   int d^3x5 <O(x6,z6) O(x5,.)> . <T3 T4 O~(x5,.)>_c  =  sum_b beta_{cb} <T3(z3) T4(z4) O(x6,z6)>_b ,      x6 = (-e)^+,  x3, x4 at the crystal.

* Left: the O-index contraction with the 2-pt structure is evaluation at z5 = a := I(x56).z6 (null), exactly as in l3.
  <T3 T4 O~>_c: conserved, 1<->2-symmetric parity-even structures at dimension 3 - Delta (l5_spin.conserved(3-Delta, J)), numerators as a
  TPS prefactor (l1_tps.expand_pref), denominators (Delta-parts in the fast exponent as in l3) in the slow exponent.
* Right: the conserved structures at dimension Delta (l5_spin.conserved(Delta, J)), the SAME basis as used for <T1 T2 O> in the functional.
* Fit over ntrip random null polarisation triples (z3, z4, z6): S_c(z) = sum_b b_{cb} r_b(z); beta = exp(logI0 - rhs_den) b.
Validation (recipe): fit residual ~ Delta^-3 (truncation) and beta independent of the crystal point nu/Delta and of the polarisations.
Usage: python3 l5_shadow_matrix.py -> cache/l5_shadow_matrix.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l1_tps as T
from l1_crystal import dot, sub, scal, cpt
import l5_spin as S


def _rnull(rng):
    y = rng.normal() + 1j*rng.normal(); u = rng.normal() + 1j*rng.normal()
    return np.array([u, y*y/u, y])


def crystal_x3x4(xi):
    a = (1 - np.sqrt(1 - xi**2))/xi
    y = np.array([1j*a, -1j*a, 0.0]); x = np.array([-1j*a, 1j*a, 0.0])
    xsq = -x[0]*x[1]
    return y, x/(-xsq)


def _setup(zv, Delta, J, x3, x4):
    w0, w1, w2 = zv[0], zv[1], zv[2]
    e = cpt(1.0, 1.0, 0.0)
    x5 = [e[0] + w0, e[1] + w1, e[2] + w2]
    x6 = cpt(-1.0, -1.0, 0.0)
    X3 = [T.const(c) for c in x3]; X4 = [T.const(c) for c in x4]
    def msq(p, q):
        r = sub(p, q); return -dot(r, r)
    d56 = msq(x5, x6)
    m35, m45, m34 = msq(X3, x5), msq(X4, x5), msq(X3, X4)
    lamF = Delta*(-T.log(d56) + 0.5*T.log(m35) + 0.5*T.log(m45))
    for i in range(3, T.NV):
        lamF = lamF - 0.5*T.mul(zv[i], zv[i])
    # O~ (dimension 3 - Delta) side denominators without the Delta-part:  X34^{-(7+Delta-J)/2} X35^{-(3+J)/2} X45^{-(3+J)/2}
    G = -((3 + J)/2)*(T.log(m35) + T.log(m45)) - ((7 + Delta - J)/2)*T.log(m34)
    x56 = sub(x5, x6)
    return lamF, G, x5, x56, X3, X4


def _avec(x56, z6):
    z6t = [T.const(c) for c in z6]
    r = T.div(dot(x56, z6t), dot(x56, x56))
    return sub(z6t, scal(2*r, x56))


def lhs_series(Delta, J, xi, trips):
    """returns logI0 (polarisation-independent) and S[t, c] (series sums, orders 0 / +1 / +2)."""
    x3, x4 = crystal_x3x4(xi)
    zs = np.zeros(T.NV, complex)
    zv = [T.var(i, 0.0) for i in range(T.NV)]
    lamF, G, *_ = _setup(zv, Delta, J, x3, x4)
    _, H = T.grad_hess(lamF)
    L = np.linalg.inv(T.csym_cholesky(-H)).T
    for _ in range(2):
        lamF2, G2, *_ = _setup(T.linsub(zs, L), Delta, J, x3, x4)
        _, H2 = T.grad_hess(lamF2)
        L = L @ np.linalg.inv(T.csym_cholesky(-H2)).T
    lamF2, G2, x5, x56, X3, X4 = _setup(T.linsub(zs, L), Delta, J, x3, x4)
    mons, solt, _ = S.conserved(3.0 - Delta, J)
    out = []
    for (z3, z4, z6) in trips:
        a = _avec(x56, z6)
        a_num = np.array([a[i][0] for i in range(3)])
        # numerators with z5 -> a (TPS, not a fixed sample): evaluate linearly: V5(a) etc. via basis expansion with TPS coefficients
        pref = _tto_at_tps_pol(X3, z3, X4, z4, x5, a, J, mons, solt)            # (nc, N)
        row = []
        for c in range(pref.shape[0]):
            S0, S2, S4 = T.expand_pref(lamF2, G2, pref[c])
            row.append((S0, S0 + S2, S0 + S2 + S4))
        out.append(row)
    logI0 = lamF[0] + G[0] + 1.5*np.log(2*np.pi) + np.log(np.linalg.det(L[:3, :3]))
    return logI0, np.array(out)                                                  # (ntrip, nc, 3)


def _tto_at_tps_pol(xa, za, xb, zb, x5, z5, J, mons, sol):
    """numerators of <T(xa,za) T(xb,zb) O(x5,z5)> with z5 a TPS vector (za, zb numeric)."""
    za_t, zb_t = S.tv(za), S.tv(zb)
    V1 = S.Vfun(xa, za_t, xb, x5); V2 = S.Vfun(xb, zb_t, x5, xa)
    V5b = [S.Vfun(x5, S.BASIS[m], xa, xb) for m in range(3)]
    H15b = [S.Hfun(xa, za_t, x5, S.BASIS[m]) for m in range(3)]
    H25b = [S.Hfun(xb, zb_t, x5, S.BASIS[m]) for m in range(3)]
    H12 = S.Hfun(xa, za_t, xb, zb_t)
    lin = lambda b: sum(S.bm(z5[m], b[m]) for m in range(3))
    V5, H15, H25 = lin(V5b), lin(H15b), lin(H25b)
    res = np.zeros((sol.shape[1], T.N), complex)
    pw = {}
    def p(nm, base, k):
        if (nm, k) not in pw:
            pw[(nm, k)] = T.const(1.0) if k == 0 else S.bm(p(nm, base, k - 1), base)
        return pw[(nm, k)]
    for i, (a1, a2, a5, b12, b15, b25) in enumerate(mons):
        v = S.bm(S.bm(S.bm(p('V1', V1, a1), p('V2', V2, a2)), S.bm(p('V5', V5, a5), p('H12', H12, b12))),
                 S.bm(p('H15', H15, b15), p('H25', H25, b25)))
        res += sol[i][:, None]*v[None]
    return res


def rhs(Delta, J, xi, trips):
    """rhs_den (log, polarisation-independent) and r[t, b] numerators of <T3 T4 O(x6, z6)> at dimension Delta, x6 = (-e)^+."""
    x3, x4 = crystal_x3x4(xi)
    x6 = np.array([-1.0, -1.0, 0.0])
    d = lambda p, q: -0.5*(p[0]*q[1] + p[1]*q[0]) + p[2]*q[2]
    m36 = d(x3 - x6, x3 - x6); m46 = d(x4 - x6, x4 - x6); m34 = -d(x3 - x4, x3 - x4)          # as l3.rhs_log ('+' flips x_{i6}^2)
    den = -((10 - Delta - J)/2)*np.log(m34 + 0j) - ((Delta + J)/2)*(np.log(m36 + 0j) + np.log(m46 + 0j))
    mons, sol, _ = S.conserved(Delta, J)
    r = []
    for (z3, z4, z6) in trips:
        v = S.tto_samples(S.tv(x3), S.tv(z3), S.tv(x4), z4[None, :], S.tv(x6), z6[None, :], J, mons, sol)[:, 0, 0, 0]
        r.append(v)
    return den, np.array(r)


def shadow_matrix(Delta, J, xi=0.5, ntrip=6, seed=7):
    rng = np.random.default_rng(seed)
    trips = [(_rnull(rng), _rnull(rng), _rnull(rng)) for _ in range(ntrip)]
    logI0, Sv = lhs_series(Delta, J, xi, trips)
    den, r = rhs(Delta, J, xi, trips)
    res = {}
    for o in range(3):
        b, *_ = np.linalg.lstsq(r, Sv[:, :, o], rcond=None)                     # r (t x b) @ b (b x c) = S (t x c)
        resid = np.linalg.norm(r @ b - Sv[:, :, o])/np.linalg.norm(Sv[:, :, o])
        res[o] = (b.T, resid)                                                     # b.T[c, b]
    return dict(logK=logI0 - den, b=res[2][0], resid=[res[o][1] for o in range(3)], b_orders=[res[o][0] for o in range(3)])


if __name__ == '__main__':
    out = []
    def pr(s=""):
        print(s, flush=True); out.append(s)
    pr("# l5_shadow_matrix: beta_{cb} = exp(logK) b_{cb} for parity-even conserved <TTO>; fit residual per Laplace order 0/1/2;")
    pr("#   stability: beta at nu/Delta = 0.4 vs 0.6 and for a second set of polarisation triples (seed)")
    for J in (0, 2, 4, 6):
        for Dl in (100.0, 300.0, 1000.0):
            A = shadow_matrix(Dl, J, 0.4, seed=7); B = shadow_matrix(Dl, J, 0.6, seed=7); Cc = shadow_matrix(Dl, J, 0.4, seed=11)
            bA = A['b']; bB = np.exp(B['logK'] - A['logK'])*B['b']; bC = np.exp(Cc['logK'] - A['logK'])*Cc['b']
            dxi = np.abs(bA - bB).max()/np.abs(bA).max(); dpol = np.abs(bA - bC).max()/np.abs(bA).max()
            # order-1 truncation vs order-2 for scale
            b1 = A['b_orders'][1]
            pr(f"  J={J} Delta={Dl:6.0f}: resid (orders 0/1/2) {A['resid'][0]:.1e}/{A['resid'][1]:.1e}/{A['resid'][2]:.1e}; "
               f"|beta(0.4)-beta(0.6)|/|beta| {dxi:.1e}; |beta(seed7)-beta(seed11)|/|beta| {dpol:.1e}; "
               f"order-2 - order-1 {np.abs(bA - b1).max()/np.abs(bA).max():.1e}")
            if J >= 4 and Dl == 1000.0:
                pr(f"     b (normalised) = {np.array2string(A['b']/np.abs(A['b']).max(), precision=6)}")
    open(os.path.join(HERE, 'cache', 'l5_shadow_matrix.txt'), 'w').write("\n".join(out) + "\n")
