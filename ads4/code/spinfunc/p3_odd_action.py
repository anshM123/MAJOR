"""
p3 (parity-odd sector, PARITY_ODD_ROUTEA.md sect. 3): the route-A' functionals C^_2 (dual q^2/V^4, f = v'-u', celestial j' = 0) and
C^_3 (dual o_H q/V^4, f = -1, celestial j' = 1) on <TTTT> blocks with the FULL conserved <TTO> structure set
  [parity-even (l5_spin.conserved, 1 or 2 for even J), parity-odd (p1_tto_odd.conserved_odd, 1 for J = 0, 2, 4, 5, 6, ...)],
as the CLS-D crystal integral expanded to relative O(1/Delta^2) (l1 engine), at either S^0 point n0 = +-1.

Integrand: identical to l6_spin_action.build_spin / l9_points.build_spin (C2, C3H), with the structure axis extended by the odd structure:
  12 side: odd numerator (eps.X_ij x even monomials, batched over the T2 / O polarisation samples) x (m12 m15 m25)^{-1/2},
  34 side (O~ at 3 - Delta): odd numerator x (m34 m35 m45)^{-1/2}            (m_ij = -x_ij^2: the conventions of l6 and of p2's LHS),
  block  M_{ab} = sum_c P_{ac} (beta^{-1})_{bc}, beta from p2_shadow_odd (block-diagonal, real).
At n0 = -1 the j' = 1 celestial sign of C^_3 is -1 (it is +1 at n0 = +1, absorbed in N^_3 as in l9).
Parity: under the transverse reflection the n0 = -1 integral equals (-1)^{#eps} times the n0 = +1 one, so the S^0 sum keeps the EE and OO
blocks (doubled) and cancels EO exactly: checked numerically by evaluating both points.
Also: a numerically checked branch diagnostic (log B vs the sum of individual logs) for the half-integer powers.
Cache: cache/p3_pts_<kind>_J<J>.pkl keyed (Delta, nu rounded 1e-9, n0); orders separate (M0, M2, M4), logC, logK, scalar (slog, sS).
"""
import os, sys, pickle, time, functools
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l1_tps as T
assert T.DEG == 6
import l1_crystal as Cr
import l5_spin as S
import l6_spin_action as A6
import p1_tto_odd as P1
import p2_shadow_odd as P2


# ------------------------------------------------------------------ epsilon by cofactors (copy of l9_points.eps_cof, l9 not imported)
def eps_cof(rows):
    """cofactor vector c_i of the 5x5 determinant with first row free: det([r0; rows]) = sum_i r0_i c_i; rows: 4 lists of 5 TPS."""
    R1, R2, R3, R4 = rows
    def m2(A, B, i, j):
        return T.mul(A[i], B[j]) - T.mul(A[j], B[i])
    P12 = {(i, j): m2(R1, R2, i, j) for i in range(5) for j in range(i + 1, 5)}
    P34 = {(i, j): m2(R3, R4, i, j) for i in range(5) for j in range(i + 1, 5)}
    c = []
    for i in range(5):
        cols = [k for k in range(5) if k != i]
        d = T.const(0.0)
        for a in range(4):
            for b in range(a + 1, 4):
                rest = [cols[k] for k in range(4) if k not in (a, b)]
                perm = [a, b] + [k for k in range(4) if k not in (a, b)]
                sgn = np.linalg.det(np.eye(4)[perm])
                d = d + sgn*T.mul(P12[(cols[a], cols[b])], P34[(rest[0], rest[1])])
        c.append((-1)**i*d)
    return c


def vec5X(x):
    return [T.const(1.0), S.ld(x, x), x[0], x[1], x[2]]


def vec5Z(x, z):
    return [T.const(0.0), 2*S.ld(x, z), z[0], z[1], z[2]]


def odd_structure_H(x0, z0, x2, x4, W2):
    """o_H(z2, z4) samples (5, 5, N) (l9_points.odd_structure, which = 'H'; doc 21.6)."""
    c = eps_cof([vec5X(x2), vec5X(x4), vec5X(x0), vec5Z(x0, z0)])
    e2b = [sum(T.mul(vec5Z(x2, S.BASIS[m])[i], c[i]) for i in range(5)) for m in range(3)]
    e4b = [-sum(T.mul(vec5Z(x4, S.BASIS[m])[i], c[i]) for i in range(5)) for m in range(3)]
    e2 = S.lin_samples(e2b, W2); e4 = S.lin_samples(e4b, W2)
    X20, X40, X24 = S.Xfun(x2, x0), S.Xfun(x4, x0), S.Xfun(x2, x4)
    n2 = T.powr(T.mul(T.mul(X24, X20), T.mul(X40, T.mul(X40, X40))), -0.5)
    n4 = T.powr(T.mul(T.mul(X24, X40), T.mul(X20, T.mul(X20, X20))), -0.5)
    H40 = S.lin_samples([S.Hfun(x4, S.BASIS[m], x0, z0) for m in range(3)], W2)
    H20 = S.lin_samples([S.Hfun(x2, S.BASIS[m], x0, z0) for m in range(3)], W2)
    return S.bm(S.bm(e2, n2)[:, None], H40[None, :]) + S.bm(H20[:, None], S.bm(e4, n4)[None, :])


# ------------------------------------------------------------------ batched parity-odd <T T O> numerators
def tto_odd_samples(xa, za, xb, zbs, x5, z5s, J, mons, sol):
    """parity-odd <T(xa,za) T(xb,zb) O(x5,z5)> NUMERATORS (p1 basis x sol; eps.X_ij included; the common denominator incl. the extra
    (X_ab X_a5 X_b5)^{1/2} NOT included), za fixed (TPS 3-vector), zb in zbs (nb x 3), z5 in z5s (n5 x 3).  Returns (1, nb, n5, N)."""
    Xa, Xb, X5 = vec5X(xa), vec5X(xb), vec5X(x5)
    Za = vec5Z(xa, za)
    Zb = [vec5Z(xb, S.BASIS[m]) for m in range(3)]
    Z5 = [vec5Z(x5, S.BASIS[k]) for k in range(3)]
    V1 = S.Vfun(xa, za, xb, x5)
    V2 = S.lin_samples([S.Vfun(xb, S.BASIS[m], x5, xa) for m in range(3)], zbs)
    V5 = S.lin_samples([S.Vfun(x5, S.BASIS[m], xa, xb) for m in range(3)], z5s)
    H12 = S.lin_samples([S.Hfun(xa, za, xb, S.BASIS[m]) for m in range(3)], zbs)
    H15 = S.lin_samples([S.Hfun(xa, za, x5, S.BASIS[m]) for m in range(3)], z5s)
    H25b = [[S.Hfun(xb, S.BASIS[m], x5, S.BASIS[k]) for k in range(3)] for m in range(3)]
    H25 = np.einsum('im,jk,mkn->ijn', zbs, z5s, np.array(H25b))
    Xab, Xa5, Xb5 = S.Xfun(xa, xb), S.Xfun(xa, x5), S.Xfun(xb, x5)
    c = eps_cof([Za, Xa, Xb, X5])                         # eps(R, Za, Xa, Xb, X5) = sum R_i c_i
    # eps(Za, Zb, Xa, Xb, X5) = -eps(Zb, Za, ...);  eps(Za, Z5, ...) = -eps(Z5, Za, ...)
    E12b = [-sum(T.mul(Zb[m][i], c[i]) for i in range(5)) for m in range(3)]
    E15b = [-sum(T.mul(Z5[k][i], c[i]) for i in range(5)) for k in range(3)]
    E12 = S.bm(S.lin_samples(E12b, zbs), Xab)                                   # (nb, N)
    E15 = S.bm(S.lin_samples(E15b, z5s), Xa5)                                   # (n5, N)
    E25m = []
    for m in range(3):
        cm = eps_cof([Zb[m], Xa, Xb, X5])                 # eps(Z5, Zb_m, Xa, Xb, X5) = sum Z5_i cm_i;  eps(Zb, Z5, ..) = -that
        E25m.append([-sum(T.mul(Z5[k][i], cm[i]) for i in range(5)) for k in range(3)])
    E25 = S.bm(np.einsum('im,jk,mkn->ijn', zbs, z5s, np.array(E25m)), Xb5)     # (nb, n5, N)
    def pw(base, kmax):
        out = [np.broadcast_to(T.const(1.0), base.shape).copy()]
        for _ in range(kmax):
            out.append(S.bm(out[-1], base))
        return out
    V1p, V2p, V5p, H12p, H15p, H25p = pw(V1, 2), pw(V2, 2), pw(V5, J), pw(H12, 2), pw(H15, 2), pw(H25, 2)
    nb, n5 = zbs.shape[0], z5s.shape[0]
    res = np.zeros((1, nb, n5, T.N), complex)
    for i, (typ, a1, a2, a5, b12, b15, b25) in enumerate(mons):
        if abs(sol[i, 0]) < 1e-14:
            continue
        A = S.bm(S.bm(V1p[a1], V2p[a2]), H12p[b12])                             # (nb, N)
        B = S.bm(V5p[a5], H15p[b15])                                             # (n5, N)
        Mv = S.bm(S.bm(A[:, None, :], B[None, :, :]), H25p[b25])                 # (nb, n5, N)
        if typ == 12:
            Mv = S.bm(Mv, E12[:, None, :])
        elif typ == 15:
            Mv = S.bm(Mv, E15[None, :, :])
        else:
            Mv = S.bm(Mv, E25)
        res[0] += sol[i, 0]*Mv
    return res


@functools.lru_cache(maxsize=None)
def struct_sets(Delta, J):
    """(even (mons, sol) or None, odd (mons, sol) or None) at Delta and 3 - Delta."""
    out = {}
    for D in (Delta, 3.0 - Delta):
        ev = S.conserved(D, J)[:2] if J % 2 == 0 else None
        mo, so, _ = P1.conserved_odd(D, J)
        od = (mo, so) if so is not None else None
        out[D] = (ev, od)
    return out


def build(zv, nu, Delta, J, kind, n0=1.0):
    lamF, G0, parts = Cr.build(zv, nu, Delta, 0, n0=n0, parts=True)
    g = parts['g']
    Pg = A6.geom(zv, n0)
    x0, x1, x2, x3, x4, e, z0 = (Pg[k] for k in ('x0', 'x1', 'x2', 'x3', 'x4', 'e', 'z0'))
    msq = A6.msq
    m12, m34, m13, m24, m14, m23 = msq(x1, x2), msq(x3, x4), msq(x1, x3), msq(x2, x4), msq(x1, x4), msq(x2, x3)
    m15, m25, m35, m45 = msq(x1, e), msq(x2, e), msq(x3, e), msq(x4, e)
    v042 = S.Vfun(x0, z0, x4, x2)
    if kind == 'C2':
        up = T.div(T.mul(m12, m34), T.mul(m13, m24)); vp = T.div(T.mul(m14, m23), T.mul(m13, m24))
        lf = T.log(vp - up)
    else:
        lf = T.const(0.0)
    Gs = (g['meas'] + g['xsq'] + g['FP'] + g['dual'] + g['cel'] + lf - 4*T.log(v042)
          - ((10 - J)/2)*T.log(m12) - (J/2)*(T.log(m15) + T.log(m25))
          - ((7 - J)/2)*T.log(m34) - ((3 + J)/2)*(T.log(m35) + T.log(m45)))
    W2, WJ = S.null_circle(2), S.null_circle(J)
    k2, kJ = S.kappa(2), S.kappa(J)
    x10, x30 = Cr.sub(x1, x0), Cr.sub(x3, x0)
    sets = struct_sets(Delta, J)
    (evD, odD), (evT, odT) = sets[Delta], sets[3.0 - Delta]
    Fs, Gts = [], []
    if evD is not None:
        Fs.append(S.tto_samples(x1, x10, x2, W2, e, WJ, J, *evD)); Gts.append(S.tto_samples(x3, x30, x4, W2, e, WJ, J, *evT))
    if odD is not None:
        ex12 = T.powr(T.mul(T.mul(m12, m15), m25), -0.5); ex34 = T.powr(T.mul(T.mul(m34, m35), m45), -0.5)
        Fs.append(S.bm(tto_odd_samples(x1, x10, x2, W2, e, WJ, J, *odD), ex12))
        Gts.append(S.bm(tto_odd_samples(x3, x30, x4, W2, e, WJ, J, *odT), ex34))
    F = np.concatenate(Fs, axis=0); Gt = np.concatenate(Gts, axis=0)
    F = S.fourier(S.fourier(F, 1, 2), 2, J); Gt = S.fourier(S.fourier(Gt, 1, 2), 2, J)
    H20 = S.lin_samples([S.Hfun(x2, S.BASIS[m], x0, z0) for m in range(3)], W2)
    H40 = S.lin_samples([S.Hfun(x4, S.BASIS[m], x0, z0) for m in range(3)], W2)
    H24b = np.array([[S.Hfun(x2, S.BASIS[m], x4, S.BASIS[k]) for k in range(3)] for m in range(3)])
    H24 = np.einsum('im,jk,mkn->ijn', W2, W2, H24b)
    inv = T.recip(T.mul(S.Xfun(x2, x0), S.Xfun(x4, x0)))
    q = S.bm(2*S.bm(H20[:, None], H40[None, :]) + S.bm(S.bm(v042, v042)[None, None], H24), inv)
    if kind == 'C2':
        D = S.bm(q, q)
    else:
        D = S.bm(q, odd_structure_H(x0, z0, x2, x4, W2))
        if n0 < 0:
            D = -D                                       # j' = 1 celestial sign at the second S^0 point
    D = S.fourier(S.fourier(D, 0, 2), 1, 2)
    na, nc = F.shape[0], Gt.shape[0]
    pref = np.zeros((na, nc, T.N), complex)
    for a in range(na):
        E = np.einsum('i,ijkn->jkn', 1/k2, S.bm(D[:, :, None, :], F[a][::-1][:, None, :, :]))
        for c in range(nc):
            Gr = Gt[c][::-1, ::-1]
            pref[a, c] = np.einsum('j,k,jkn->n', 1/k2, 1/kJ, S.bm(E, Gr))
    diag = dict(logB=parts['lB'][0], logm=[T.log(mm)[0] for mm in (m35, m45, m12, m15, m25, m34)])
    return lamF, Gs, pref, diag


def action(nu, Delta, J, kind='C2', n0=1.0):
    zs, a = Cr.saddle(nu, Delta)
    zv = [T.var(i, zs[i]) for i in range(T.NV)]
    lamF = Cr.build(zv, nu, Delta, 0, n0=n0)[0]
    _, HF = T.grad_hess(lamF)
    L = np.linalg.inv(T.csym_cholesky(-HF)).T
    for _ in range(2):
        lamF2 = Cr.build(T.linsub(zs, L), nu, Delta, 0, n0=n0)[0]
        _, H2 = T.grad_hess(lamF2)
        L = L @ np.linalg.inv(T.csym_cholesky(-H2)).T
    lamF2, G2, pref, diag = build(T.linsub(zs, L), nu, Delta, J, kind, n0)
    na, nc = pref.shape[:2]
    Ser = np.zeros((na, nc, 3), complex)
    for i in range(na):
        for c in range(nc):
            Ser[i, c] = T.expand_pref(lamF2, G2, pref[i, c])
    b, logK, resid = P2.shadow_odd(Delta, J)
    M = np.einsum('aco,bc->abo', Ser, np.linalg.inv(b))
    logC = lamF[0] + G2[0] + np.log(np.linalg.det(L))
    br = diag['logB'] - (diag['logm'][0] + diag['logm'][1] + diag['logm'][2] - diag['logm'][3] - diag['logm'][4] - diag['logm'][5])
    return dict(M0=M[:, :, 0], M2=M[:, :, 1], M4=M[:, :, 2], logC=logC, logK=logK, shadow_resid=resid, branch=br)


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def _fn(kind, J):
    return os.path.join(HERE, 'cache', 'p3_pts_%s_J%d.pkl' % (kind, J))


def load(kind, J):
    fn = _fn(kind, J)
    return pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}


def get(kind, J, Delta, nus, n0=1.0, verbose=True):
    res = load(kind, J)
    out = []
    for nu in nus:
        key = (float(Delta), round(float(nu), 9), float(n0))
        if key not in res:
            t = time.time()
            r = action(float(nu), float(Delta), J, kind, n0)
            s = Cr.laplace_pref(float(nu), float(Delta), J, n0=n0)
            r.update(slog=s['logY0'] - 3.5*np.log(2*np.pi), sS=np.array(s['S'][:3]), m=mphys(Delta, J))
            res = load(kind, J); res[key] = r
            tmp = _fn(kind, J) + '.tmp%d' % os.getpid()
            pickle.dump(res, open(tmp, 'wb')); os.replace(tmp, _fn(kind, J))
            if verbose:
                print(f"{kind} J={J} Delta={Delta} nu={nu:.4f} n0={n0:+.0f}: {time.time() - t:.1f}s", flush=True)
        out.append(res[key])
    return out


if __name__ == '__main__':
    # usage: python3 p3_odd_action.py kind J Delta n0 x1 [x2 ...]   (x = nu/m)
    kind, J, Dl, n0 = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    m = mphys(Dl, J)
    get(kind, J, Dl, [float(x)*m for x in sys.argv[5:]], n0=n0)
