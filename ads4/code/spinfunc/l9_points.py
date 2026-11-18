"""
l9_points (M6 round 5, doc 23): graviton crystal actions with the NEXT Laplace order K3 (eps^6; l9_tps, DEG = 8), for
  kind 'C2' : C^_2  (dual q^2/V^4, f = v'-u', celestial j' = 0)            -- as l6_spin_action.build_spin
  kind 'C3H': C^_3  (dual o_H q/V^4, f = -1, celestial j' = 1 = j' = 0 x locally constant sign)
  kind 'C3V': C^_3  with the V-type completion o_V of the odd structure (same crystal image; doc 21.6 design freedom)
     o_H = e(Z2,X2,X4,X0,Z0) H40/sqrt(X24 X20 X40^3) + e(Z4,X4,X2,X0,Z0) H20/sqrt(X24 X40 X20^3)                      (g6, doc 21.6)
     o_V = V042 [ e(Z2,X2,X4,X0,Z0) V4,02/sqrt(X24 X20 X40^3) - V2,04 e(Z4,X4,X2,X0,Z0)/sqrt(X24 X40 X20^3) ]         (homogeneity as o_H)
  Both n0 = +-1 give the same value for every kind (o and the j' = 1 sign are both odd under the S^0 reflection); n0 = +1 is used.
Orders kept separate: M0, M2, M4, M6 (2x2 / 1x1 in the conserved structure basis of l5_spin.conserved), logC, logK; scalar k = 2 crystal
(prefactor form) S0..S6.  The shadow matrix is computed at K3 as well (l1_tps.expand_pref patched to the exact Hermite engine, third
entry = S4 + S6), cached per (Delta, J) in cache/l9_shadow.pkl.
Cache: cache/l9_pts_<kind>_J<J>.pkl keyed (Delta, nu rounded 1e-9).  Resumable; several processes may write (re-load before dump).
Usage: L1_DEG=8 python3 l9_points.py kind J Delta nu1 [nu2 ...]
"""
import os, sys, pickle, time
os.environ.setdefault('L1_DEG', '8')                  # set L1_DEG=6 for the K2-only (fast) variant
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l1_tps as T
K3 = T.DEG >= 8                                   # L1_DEG=6: K2 only (fast; M6 = 0), separate caches (suffix _k2)
if K3:
    import l9_tps as K9
    T.expand_pref = lambda lamF, G, q=None: (lambda s: (s[0], s[1], s[2] + s[3]))(K9.expand_pref6(lamF, G, q))   # shadow at K3
    EXP6 = K9.expand_pref6
else:
    EXP6 = lambda lamF, G, q=None: tuple(T.expand_pref(lamF, G, q)) + (0.0,)
SFX = '' if K3 else '_k2'
import l1_crystal as Cr
import l5_spin as S
import l5_shadow_matrix as SM
import l6_spin_action as A6

KINDS = ('C2', 'C3H', 'C3V')


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def shadow(Delta, J):
    fn = os.path.join(HERE, 'cache', 'l9_shadow%s.pkl' % SFX)
    db = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    key = (float(Delta), int(J))
    if key not in db:
        r = SM.shadow_matrix(Delta, J, 0.5, ntrip=6, seed=7)
        db = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
        db[key] = (r['b'], r['logK'], r['resid'])
        pickle.dump(db, open(fn + '.tmp%d' % os.getpid(), 'wb')); os.replace(fn + '.tmp%d' % os.getpid(), fn)
    return db[key]


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


def odd_structure(x0, z0, x2, x4, W2, which):
    """o(z2, z4) samples on the null circles (5, 5, N) for which = 'H' or 'V'."""
    c = eps_cof([vec5X(x2), vec5X(x4), vec5X(x0), vec5Z(x0, z0)])              # eps(Z2, X2, X4, X0, Z0) = sum Z2_i c_i
    e2b = [sum(T.mul(vec5Z(x2, S.BASIS[m])[i], c[i]) for i in range(5)) for m in range(3)]
    e4b = [-sum(T.mul(vec5Z(x4, S.BASIS[m])[i], c[i]) for i in range(5)) for m in range(3)]  # eps(Z4,X4,X2,X0,Z0) = -eps(Z4,X2,X4,..)
    e2 = S.lin_samples(e2b, W2); e4 = S.lin_samples(e4b, W2)                      # (5, N)
    X20, X40, X24 = S.Xfun(x2, x0), S.Xfun(x4, x0), S.Xfun(x2, x4)
    n2 = T.powr(T.mul(T.mul(X24, X20), T.mul(X40, T.mul(X40, X40))), -0.5)
    n4 = T.powr(T.mul(T.mul(X24, X40), T.mul(X20, T.mul(X20, X20))), -0.5)
    if which == 'H':
        H40 = S.lin_samples([S.Hfun(x4, S.BASIS[m], x0, z0) for m in range(3)], W2)
        H20 = S.lin_samples([S.Hfun(x2, S.BASIS[m], x0, z0) for m in range(3)], W2)
        o = S.bm(S.bm(e2, n2)[:, None], H40[None, :]) + S.bm(H20[:, None], S.bm(e4, n4)[None, :])
    else:
        v042 = S.Vfun(x0, z0, x4, x2)
        V402 = S.lin_samples([S.Vfun(x4, S.BASIS[m], x0, x2) for m in range(3)], W2)
        V204 = S.lin_samples([S.Vfun(x2, S.BASIS[m], x0, x4) for m in range(3)], W2)
        o = S.bm(v042, S.bm(S.bm(e2, n2)[:, None], V402[None, :]) - S.bm(V204[:, None], S.bm(e4, n4)[None, :]))
    return o                                                                      # (5[z2], 5[z4], N)


def build_spin(zv, nu, Delta, J, mons, sol, solt, kind):
    lamF, G0, parts = Cr.build(zv, nu, Delta, 0, n0=1.0, parts=True)
    g = parts['g']
    P = A6.geom(zv, 1.0)
    x0, x1, x2, x3, x4, e, z0 = (P[k] for k in ('x0', 'x1', 'x2', 'x3', 'x4', 'e', 'z0'))
    msq = A6.msq
    m12, m34, m13, m24, m14, m23 = msq(x1, x2), msq(x3, x4), msq(x1, x3), msq(x2, x4), msq(x1, x4), msq(x2, x3)
    m15, m25, m35, m45 = msq(x1, e), msq(x2, e), msq(x3, e), msq(x4, e)
    v042 = S.Vfun(x0, z0, x4, x2)
    if kind == 'C2':
        up = T.div(T.mul(m12, m34), T.mul(m13, m24)); vp = T.div(T.mul(m14, m23), T.mul(m13, m24))
        lf = T.log(vp - up)
    else:
        lf = T.const(0.0)                                                         # f = -1 (constant sign: absorbed in N^_3)
    Gs = (g['meas'] + g['xsq'] + g['FP'] + g['dual'] + g['cel'] + lf - 4*T.log(v042)
          - ((10 - J)/2)*T.log(m12) - (J/2)*(T.log(m15) + T.log(m25))
          - ((7 - J)/2)*T.log(m34) - ((3 + J)/2)*(T.log(m35) + T.log(m45)))
    W2, WJ = S.null_circle(2), S.null_circle(J)
    k2, kJ = S.kappa(2), S.kappa(J)
    x10, x30 = Cr.sub(x1, x0), Cr.sub(x3, x0)
    F = S.tto_samples(x1, x10, x2, W2, e, WJ, J, mons, sol)
    Gt = S.tto_samples(x3, x30, x4, W2, e, WJ, J, mons, solt)
    F = S.fourier(S.fourier(F, 1, 2), 2, J); Gt = S.fourier(S.fourier(Gt, 1, 2), 2, J)
    H20 = S.lin_samples([S.Hfun(x2, S.BASIS[m], x0, z0) for m in range(3)], W2)
    H40 = S.lin_samples([S.Hfun(x4, S.BASIS[m], x0, z0) for m in range(3)], W2)
    H24b = np.array([[S.Hfun(x2, S.BASIS[m], x4, S.BASIS[k]) for k in range(3)] for m in range(3)])
    H24 = np.einsum('im,jk,mkn->ijn', W2, W2, H24b)
    inv = T.recip(T.mul(S.Xfun(x2, x0), S.Xfun(x4, x0)))
    q = S.bm(2*S.bm(H20[:, None], H40[None, :]) + S.bm(S.bm(v042, v042)[None, None], H24), inv)
    D = S.bm(q, q) if kind == 'C2' else S.bm(q, odd_structure(x0, z0, x2, x4, W2, kind[-1]))
    D = S.fourier(S.fourier(D, 0, 2), 1, 2)
    na, nc = F.shape[0], Gt.shape[0]
    pref = np.zeros((na, nc, T.N), complex)
    for a in range(na):
        E = np.einsum('i,ijkn->jkn', 1/k2, S.bm(D[:, :, None, :], F[a][::-1][:, None, :, :]))
        for c in range(nc):
            Gr = Gt[c][::-1, ::-1]
            pref[a, c] = np.einsum('j,k,jkn->n', 1/k2, 1/kJ, S.bm(E, Gr))
    return lamF, Gs, pref


def action(nu, Delta, J, kind='C2'):
    mons, sol, _ = S.conserved(Delta, J)
    _, solt, _ = S.conserved(3.0 - Delta, J)
    zs, a = Cr.saddle(nu, Delta)
    zv = [T.var(i, zs[i]) for i in range(T.NV)]
    lamF = Cr.build(zv, nu, Delta, 0)[0]
    _, HF = T.grad_hess(lamF)
    L = np.linalg.inv(T.csym_cholesky(-HF)).T
    for _ in range(2):
        lamF2 = Cr.build(T.linsub(zs, L), nu, Delta, 0)[0]
        _, H2 = T.grad_hess(lamF2)
        L = L @ np.linalg.inv(T.csym_cholesky(-H2)).T
    lamF2, G2, pref = build_spin(T.linsub(zs, L), nu, Delta, J, mons, sol, solt, kind)
    na, nc = pref.shape[:2]
    Ser = np.zeros((na, nc, 4), complex)
    for i in range(na):
        for c in range(nc):
            Ser[i, c] = EXP6(lamF2, G2, pref[i, c])
    b, logK, resid = shadow(Delta, J)
    M = np.einsum('aco,bc->abo', Ser, np.linalg.inv(b))
    logC = lamF[0] + G2[0] + np.log(np.linalg.det(L))
    return dict(M=M, logC=logC, logK=logK, shadow_resid=resid)


def scalar6(nu, Delta, J):
    zs, a = Cr.saddle(nu, Delta)
    zv = [T.var(i, zs[i]) for i in range(T.NV)]
    lamF, G, Q = Cr.build(zv, nu, Delta, J, pref=True)
    _, HF = T.grad_hess(lamF)
    L = np.linalg.inv(T.csym_cholesky(-HF)).T
    for _ in range(2):
        lamF2, G2, Q2 = Cr.build(T.linsub(zs, L), nu, Delta, J, pref=True)
        _, H2 = T.grad_hess(lamF2)
        L = L @ np.linalg.inv(T.csym_cholesky(-H2)).T
    lamF2, G2, Q2 = Cr.build(T.linsub(zs, L), nu, Delta, J, pref=True)
    S6 = EXP6(lamF2, G2, Q2)
    logY0 = lamF[0] + G[0] + 3.5*np.log(2*np.pi) + np.log(np.linalg.det(L)) - Cr.log_beta(Delta, J)
    return logY0 - 3.5*np.log(2*np.pi), np.array(S6)


def _fn(kind, J):
    return os.path.join(HERE, 'cache', 'l9_pts_%s_J%d%s.pkl' % (kind, J, SFX))


def load(kind, J):
    fn = _fn(kind, J)
    return pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}


def get(kind, J, Delta, nus, verbose=True):
    res = load(kind, J)
    if kind == 'C2' and not K3:                                     # K2: reuse the l8/l7 C2 records (same engine, same orders)
        import l8_analysis as L8
        R8 = L8.records(J); D8 = {}
        for k in R8:
            D8.setdefault(k[0], []).append(k[1])
        for nu in nus:                                              # match nu to 1e-7 relative (l7 keys use l7's own m)
            key = (float(Delta), round(float(nu), 9))
            if key not in res and float(Delta) in D8:
                c = min(D8[float(Delta)], key=lambda v: abs(v - nu))
                if abs(c - nu) < 1e-7*nu:
                    r = R8[(float(Delta), c)]; res[key] = dict(r, M6=0*r['M0'], sS=np.append(r['sS'], 0.0))
    out = []
    for nu in nus:
        key = (float(Delta), round(float(nu), 9))
        if key not in res:
            t = time.time()
            r = action(float(nu), float(Delta), J, kind)
            sl, sS = scalar6(float(nu), float(Delta), J)
            M = r['M']
            rec = dict(M0=M[:, :, 0], M2=M[:, :, 1], M4=M[:, :, 2], M6=M[:, :, 3], logC=r['logC'], logK=r['logK'],
                       slog=sl, sS=sS, m=mphys(Delta, J), shadow_resid=r['shadow_resid'])
            res = load(kind, J); res[key] = rec
            tmp = _fn(kind, J) + '.tmp%d' % os.getpid()
            pickle.dump(res, open(tmp, 'wb')); os.replace(tmp, _fn(kind, J))
            if verbose:
                print(f"{kind} J={J} Delta={Delta} nu={nu:.4f}: {time.time() - t:.1f}s", flush=True)
        out.append(res[key])
    return out


if __name__ == '__main__':
    kind = sys.argv[1]; J = int(sys.argv[2]); Dl = float(sys.argv[3])
    get(kind, J, Dl, [float(v) for v in sys.argv[4:]])
