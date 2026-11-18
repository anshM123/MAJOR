"""
l6 (M6 step 2, items 2-3 of the 20.9 recipe): the route-A' functional C^_2 (j' = 0, (2,4) structure q^2, subtraction f+_{2;2} = v'-u') on a
parity-even <TTTT> s-channel block G_{(a'b')}[Delta, J], as the CLS-D crystal integral expanded to relative O(1/Delta^2) with the l1 engine.

Integrand (CLS-D eq. spinningfunctional_definition / Psi^+, gauge x3 = -x1 = y, x4 = -x2 = x/(-x^2), x5 = e, d = 3, delta functions solved
as in l1_crystal):
   [measure, (-x^2)^{-2d}, FP]                                                  (l1_crystal, unchanged)
 x [scalar dual, Delta~_{2,4} = 0, J = -1] x q(z2,z4)^2 / V_{0,42}^4            (q of g5_covariant_Q / doc 17; the V^-4 restores the
                                                                                  homogeneity of the J = 3 dual, see doc 21.2)
 x [celestial j' = 0, delta_i = 2]                                              (l1_crystal, unchanged: light-transform weights cancel)
 x (v' - u')                                                                     (f+_{2;2})
 x <T1(z1 = x10) T2 O(e)>_{a'} . beta^{-1} . <O~(e) T3(z3 = x30) T4>_c           (l5_spin numerators; denominators: Delta-parts in the
                                                                                  CLS-D fast exponent, rest in the slow exponent)
 with the T2, T4 and O indices contracted (l5_spin Fourier contraction on null circles).
Everything polarisation-dependent is one prefactor P_{a'c} (TPS), expanded with l1_tps.expand_pref; M_{a'b'} = sum_c P-series_{a'c} (beta^-1)_{b'c}.
Output of action(): logC (common, complex), Ser[a', c, order], b, logK (shadow matrix), M[order] (a' x b') up to exp(logC - logK),
and the scalar reference (l1_crystal.laplace_pref at the same nu, Delta, J) for ratios.
"""
import os, sys, functools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l1_tps as T
import l1_crystal as Cr
import l5_spin as S
import l5_shadow_matrix as SM


def geom(zv, n0=1.0):
    """gauge-fixed positions as TPS 3-vectors (copy of l1_crystal.build's geometry)."""
    xt0, ux, vx, uy, vy, xp, yp = zv
    x = [ux, vx, xp]; y = [uy, vy, yp]
    uv = T.mul(uy, vy); yp2 = T.mul(yp, yp); one = T.const(1.0)
    y0sq = T.div(T.mul(uv, T.mul(uv, one - T.mul(xt0, xt0)) - yp2), uv - yp2)
    y0 = n0*T.powr(y0sq, 0.5)
    u0 = T.mul(uy, xt0) + T.div(T.mul(y0, yp), vy)
    v0 = -T.mul(vy, xt0) + T.div(T.mul(y0, yp), uy)
    x0 = [u0, v0, y0]
    xsq = Cr.dot(x, x)
    x4 = Cr.scal(T.recip(-xsq), x); x2 = Cr.neg(x4); x1 = Cr.neg(y); x3 = y
    return dict(x0=x0, x1=x1, x2=x2, x3=x3, x4=x4, e=Cr.cpt(1.0, 1.0, 0.0), z0=Cr.cpt(1.0, 0.0, 0.0))


def msq(p, q):
    r = Cr.sub(p, q); return -Cr.dot(r, r)


@functools.lru_cache(maxsize=None)
def shadow(Delta, J):
    r = SM.shadow_matrix(Delta, J, 0.5, ntrip=6, seed=7)
    return r['b'], r['logK'], r['resid'][2]


def build_spin(zv, nu, Delta, J, mons, sol, solt, n0=1.0):
    lamF, G0, parts = Cr.build(zv, nu, Delta, 0, n0=n0, parts=True)
    g = parts['g']
    P = geom(zv, n0)
    x0, x1, x2, x3, x4, e, z0 = (P[k] for k in ('x0', 'x1', 'x2', 'x3', 'x4', 'e', 'z0'))
    m12, m34, m13, m24, m14, m23 = msq(x1, x2), msq(x3, x4), msq(x1, x3), msq(x2, x4), msq(x1, x4), msq(x2, x3)
    m15, m25, m35, m45 = msq(x1, e), msq(x2, e), msq(x3, e), msq(x4, e)
    up = T.div(T.mul(m12, m34), T.mul(m13, m24)); vp = T.div(T.mul(m14, m23), T.mul(m13, m24))
    v042 = S.Vfun(x0, z0, x4, x2)
    Gs = (g['meas'] + g['xsq'] + g['FP'] + g['dual'] + g['cel'] + T.log(vp - up) - 4*T.log(v042)
          - ((10 - J)/2)*T.log(m12) - (J/2)*(T.log(m15) + T.log(m25))
          - ((7 - J)/2)*T.log(m34) - ((3 + J)/2)*(T.log(m35) + T.log(m45)))
    # polarisation sums
    W2, WJ = S.null_circle(2), S.null_circle(J)
    k2, kJ = S.kappa(2), S.kappa(J)
    x10, x30 = Cr.sub(x1, x0), Cr.sub(x3, x0)
    F = S.tto_samples(x1, x10, x2, W2, e, WJ, J, mons, sol)                  # (na, 5, n5, N)
    Gt = S.tto_samples(x3, x30, x4, W2, e, WJ, J, mons, solt)                # (nc, 5, n5, N)
    F = S.fourier(S.fourier(F, 1, 2), 2, J); Gt = S.fourier(S.fourier(Gt, 1, 2), 2, J)
    H20 = S.lin_samples([S.Hfun(x2, S.BASIS[m], x0, z0) for m in range(3)], W2)     # (5, N)
    H40 = S.lin_samples([S.Hfun(x4, S.BASIS[m], x0, z0) for m in range(3)], W2)
    H24b = np.array([[S.Hfun(x2, S.BASIS[m], x4, S.BASIS[k]) for k in range(3)] for m in range(3)])
    H24 = np.einsum('im,jk,mkn->ijn', W2, W2, H24b)                                # (5, 5, N)
    inv = T.recip(T.mul(S.Xfun(x2, x0), S.Xfun(x4, x0)))
    q = S.bm(2*S.bm(H20[:, None], H40[None, :]) + S.bm(S.bm(v042, v042)[None, None], H24), inv)
    D = S.bm(q, q)                                                                  # (5, 5, N)
    D = S.fourier(S.fourier(D, 0, 2), 1, 2)
    na, nc = F.shape[0], Gt.shape[0]
    pref = np.zeros((na, nc, T.N), complex)
    for a in range(na):
        # E[n4, n5] = sum_n2 D[n2, n4] F_a[-n2, n5] / k2[n2]
        E = np.einsum('i,ijkn->jkn', 1/k2, S.bm(D[:, :, None, :], F[a][::-1][:, None, :, :]))
        for c in range(nc):
            Gr = Gt[c][::-1, ::-1]                                                  # Gt_c[-n4, -n5]
            pref[a, c] = np.einsum('j,k,jkn->n', 1/k2, 1/kJ, S.bm(E, Gr))
    return lamF, Gs, pref


def action(nu, Delta, J, with_scalar=True, n0=1.0):
    mons, sol, _ = S.conserved(Delta, J)
    _, solt, _ = S.conserved(3.0 - Delta, J)
    zs, a = Cr.saddle(nu, Delta)
    zv = [T.var(i, zs[i]) for i in range(T.NV)]
    lamF = Cr.build(zv, nu, Delta, 0, n0=n0)[0]
    _, HF = T.grad_hess(lamF)
    L = np.linalg.inv(T.csym_cholesky(-HF)).T
    for _ in range(2):
        lamF2 = Cr.build(T.linsub(zs, L), nu, Delta, 0, n0=n0)[0]
        _, H2 = T.grad_hess(lamF2)
        L = L @ np.linalg.inv(T.csym_cholesky(-H2)).T
    lamF2, G2, pref = build_spin(T.linsub(zs, L), nu, Delta, J, mons, sol, solt, n0)
    na, nc = pref.shape[:2]
    Ser = np.zeros((na, nc, 3), complex)
    for i in range(na):
        for c in range(nc):
            S0, S2, S4 = T.expand_pref(lamF2, G2, pref[i, c])
            Ser[i, c] = (S0, S0 + S2, S0 + S2 + S4)
    b, logK, resid = shadow(Delta, J)
    binv = np.linalg.inv(b)                                                        # (b' x c)
    M = np.einsum('aco,bc->abo', Ser, binv)                                        # (a', b', order)
    logC = lamF[0] + G2[0] + np.log(np.linalg.det(L))
    out = dict(logC=logC, logK=logK, Ser=Ser, b=b, M=M, a=a, shadow_resid=resid)
    if with_scalar:
        r = Cr.laplace_pref(nu, Delta, J)
        out['scal_log'] = r['logY0'] - 3.5*np.log(2*np.pi)
        out['scal_S'] = np.array([r['S'][0], r['S'][0] + r['S'][1], r['S'][0] + r['S'][1] + r['S'][2]])
    return out
