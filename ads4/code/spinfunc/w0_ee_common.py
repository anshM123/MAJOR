"""
w0 (EE_NEARGAP.md): the EE graviton integrand (C^_EE = Psi_2[EE], j' = 0, f = v'-u', celestial j' = 0) as module-parametrised factories,
so that the SAME code runs on the TPS engine (l1_tps; K2 Laplace rows) and on the point-wise backend (n1_pw; direct RQMC):
    q = (2 H20 H40 + V042^2 H24)/(X20 X40)      ->      E = V042^2 (H24 - 2 V_{2,04} V_{4,02})/(X20 X40),   D = E (x) E
(u1/EE_FUNCTIONAL §1: E = q at the crystal point).  Even J: the l9_points.build_spin integrand (as u5.build_spin_EE, which it
reproduces); odd J: the p3_odd_action.build integrand (full even + odd conserved <TTO> set), with E in place of q in the C2 line.
Project modules are imported read-only.
"""
import numpy as np


def e_factor(T, S, x0, z0, x2, x4, W2):
    """E(z2, z4) samples (5, 5, N) on the spin-2 null circles."""
    v042 = S.Vfun(x0, z0, x4, x2)
    H24b = np.array([[S.Hfun(x2, S.BASIS[m], x4, S.BASIS[k]) for k in range(3)] for m in range(3)])
    H24 = np.einsum('im,jk,mkn->ijn', W2, W2, H24b)
    V204 = S.lin_samples([S.Vfun(x2, S.BASIS[m], x0, x4) for m in range(3)], W2)
    V402 = S.lin_samples([S.Vfun(x4, S.BASIS[m], x0, x2) for m in range(3)], W2)
    inv = T.recip(T.mul(S.Xfun(x2, x0), S.Xfun(x4, x0)))
    return S.bm(S.bm(S.bm(v042, v042)[None, None], H24 - 2*S.bm(V204[:, None], V402[None, :])), inv)


def _gs(T, g, x0, x1, x2, x3, x4, e, z0, J, A6, S, c2=True):
    msq = A6.msq
    m12, m34, m13, m24, m14, m23 = msq(x1, x2), msq(x3, x4), msq(x1, x3), msq(x2, x4), msq(x1, x4), msq(x2, x3)
    m15, m25, m35, m45 = msq(x1, e), msq(x2, e), msq(x3, e), msq(x4, e)
    v042 = S.Vfun(x0, z0, x4, x2)
    up = T.div(T.mul(m12, m34), T.mul(m13, m24)); vp = T.div(T.mul(m14, m23), T.mul(m13, m24))
    lf = T.log(vp - up)
    Gs = (g['meas'] + g['xsq'] + g['FP'] + g['dual'] + g['cel'] + lf - 4*T.log(v042)
          - ((10 - J)/2)*T.log(m12) - (J/2)*(T.log(m15) + T.log(m25))
          - ((7 - J)/2)*T.log(m34) - ((3 + J)/2)*(T.log(m35) + T.log(m45)))
    return Gs, (m12, m34, m15, m25, m35, m45)


def _contract(S, T, D, F, Gt, J):
    k2, kJ = S.kappa(2), S.kappa(J)
    D = S.fourier(S.fourier(D, 0, 2), 1, 2)
    na, nc = F.shape[0], Gt.shape[0]
    pref = np.zeros((na, nc, T.N), complex)
    for a in range(na):
        Em = np.einsum('i,ijkn->jkn', 1/k2, S.bm(D[:, :, None, :], F[a][::-1][:, None, :, :]))
        for c in range(nc):
            Gr = Gt[c][::-1, ::-1]
            pref[a, c] = np.einsum('j,k,jkn->n', 1/k2, 1/kJ, S.bm(Em, Gr))
    return pref


def make_build_even_EE(T, S, Cr, A6):
    """build_spin(zv, nu, Delta, J, mons, sol, solt, kind) -> (lamF, Gs, pref) for kind 'EE' (l9_points signature)."""
    def build(zv, nu, Delta, J, mons, sol, solt, kind='EE'):
        lamF, G0, parts = Cr.build(zv, nu, Delta, 0, n0=1.0, parts=True)
        Pg = A6.geom(zv, 1.0)
        x0, x1, x2, x3, x4, e, z0 = (Pg[k] for k in ('x0', 'x1', 'x2', 'x3', 'x4', 'e', 'z0'))
        Gs, _ = _gs(T, parts['g'], x0, x1, x2, x3, x4, e, z0, J, A6, S)
        W2, WJ = S.null_circle(2), S.null_circle(J)
        x10, x30 = Cr.sub(x1, x0), Cr.sub(x3, x0)
        F = S.tto_samples(x1, x10, x2, W2, e, WJ, J, mons, sol)
        Gt = S.tto_samples(x3, x30, x4, W2, e, WJ, J, mons, solt)
        F = S.fourier(S.fourier(F, 1, 2), 2, J); Gt = S.fourier(S.fourier(Gt, 1, 2), 2, J)
        Eb = e_factor(T, S, x0, z0, x2, x4, W2)
        return lamF, Gs, _contract(S, T, S.bm(Eb, Eb), F, Gt, J)
    return build


def make_build_odd_EE(P3):
    """p3_odd_action.build replacement: kind 'EE' -> D = E (x) E with the C2 (f = v'-u') line; other kinds -> the original."""
    T, S, Cr, A6 = P3.T, P3.S, P3.Cr, P3.A6
    orig = P3.build

    def build(zv, nu, Delta, J, kind, n0=1.0):
        if kind != 'EE':
            return orig(zv, nu, Delta, J, kind, n0)
        lamF, G0, parts = Cr.build(zv, nu, Delta, 0, n0=n0, parts=True)
        Pg = A6.geom(zv, n0)
        x0, x1, x2, x3, x4, e, z0 = (Pg[k] for k in ('x0', 'x1', 'x2', 'x3', 'x4', 'e', 'z0'))
        Gs, (m12, m34, m15, m25, m35, m45) = _gs(T, parts['g'], x0, x1, x2, x3, x4, e, z0, J, A6, S)
        W2, WJ = S.null_circle(2), S.null_circle(J)
        x10, x30 = Cr.sub(x1, x0), Cr.sub(x3, x0)
        sets = P3.struct_sets(Delta, J)
        (evD, odD), (evT, odT) = sets[Delta], sets[3.0 - Delta]
        Fs, Gts = [], []
        if evD is not None:
            Fs.append(S.tto_samples(x1, x10, x2, W2, e, WJ, J, *evD)); Gts.append(S.tto_samples(x3, x30, x4, W2, e, WJ, J, *evT))
        if odD is not None:
            ex12 = T.powr(T.mul(T.mul(m12, m15), m25), -0.5); ex34 = T.powr(T.mul(T.mul(m34, m35), m45), -0.5)
            Fs.append(S.bm(P3.tto_odd_samples(x1, x10, x2, W2, e, WJ, J, *odD), ex12))
            Gts.append(S.bm(P3.tto_odd_samples(x3, x30, x4, W2, e, WJ, J, *odT), ex34))
        F = np.concatenate(Fs, axis=0); Gt = np.concatenate(Gts, axis=0)
        F = S.fourier(S.fourier(F, 1, 2), 2, J); Gt = S.fourier(S.fourier(Gt, 1, 2), 2, J)
        Eb = e_factor(T, S, x0, z0, x2, x4, W2)
        pref = _contract(S, T, S.bm(Eb, Eb), F, Gt, J)
        diag = dict(logB=parts['lB'][0], logm=[T.log(mm)[0] for mm in (m35, m45, m12, m15, m25, m34)])
        return lamF, Gs, pref, diag
    return build
