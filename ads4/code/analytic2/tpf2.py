"""
tpf2.py: Lemma F with first-order (affine) frames.  ANALYTIC_SLIVER2 §3.

Same setting as tpframes.py (b = 0, s = 1; box sigma in [s_lo, s_hi], tau in [t_lo, t_hi], a in [a_lo, a_hi], iota = 1/N in
[i_lo, i_hi]; last segment k = k_A .. n-1 with m = N - k = a + K - i, i = 0..K-1).  Parameters p = (tau, a, iota), centre p0.
Frames are AFFINE in p:  F_i(p) = F_i^0 + sum_l F_i^l (p_l - p0_l)  (point matrices from the variational equations at p0, so that
F_{i+1}(p) ~ M_i(p) F_i(p) to second order).  With Q_i(p) = F_i(p)^{-1}, c'_i = Q_i(p) v_i(p):
    c'_{i+1} = E_i(p) c'_i,   E_i(p) = Q_{i+1}(p) M_i(p) F_i(p),
and by the mean value theorem, for p in the box B:
    E_i(p) in E_i(p0) + sum_l dE_i/dp_l (B) [-D_l, D_l],
    dE/dp_l = -Q_{i+1} F_{i+1}^l Q_{i+1} M F_i + Q_{i+1} (dM/dp_l) F_i + Q_{i+1} M F_i^l   (all evaluated as Arb balls over B).
sigma enters only through alpha = sigma + i tau inside M (zeroth order: sigma is a ball, sigma-width charged directly).
Then  |S^_n| <= Phi_I ||Q_0 P_{k_A}|| prod_i [ (d_{k+1}/d_k)^{s_lo} ||E_i||_inf ] ||row_1 F_K||_1   (Phi_I: tpframes.phase1_log).
"""
import math, sys
import numpy as np
sys.path.insert(0, '.')
from flint import arb, acb, ctx, acb_mat
from tpframes import phase1_log, up, arb_iv, _inf_norm
ctx.prec = 96


def _eps(m, io):
    return 2*(2 - (2*m - 1)*io)/(m*(2 - (m - 1)*io))


def _deps_dio(m, io):
    B = 2 - (m - 1)*io
    return -4/(B*B)


def _deps_dm(m, io):
    A = 2 - (2*m - 1)*io; B = 2 - (m - 1)*io
    return 2*(-2*io*m*B - A*(B - m*io))/(m*B)**2


def _M(w):
    return np.array([[2*w, 1], [1, 0]], dtype=complex)


def _dM(dw):
    return np.array([[2*dw, 0], [0, 0]], dtype=complex)


def build_frames(sig, tau, a, io, K, i_s):
    """F0[i], F1[i][l] (l = 0: tau, 1: a, 2: iota) for i = 0..K, plus P0 = P_{k_A}(p0)."""
    al = complex(sig, tau)
    Ms = []; dMs = []
    for i in range(K):
        m = a + K - i
        e = _eps(m, io)
        w = e*al/2
        Ms.append(_M(w))
        dMs.append([_dM(e*1j/2), _dM(_deps_dm(m, io)*al/2), _dM(_deps_dio(m, io)*al/2)])
    w0 = _eps(a + K, io)*al/2; r0 = np.sqrt(1 + w0*w0)
    P0 = np.array([[w0 + r0, w0 - r0], [1, 1]], dtype=complex)
    Phi = [P0.copy()]; dPhi = [[np.zeros((2, 2), complex) for _ in range(3)]]
    for i in range(K):
        Phi.append(Ms[i] @ Phi[i])
        dPhi.append([dMs[i][l] @ Phi[i] + Ms[i] @ dPhi[i][l] for l in range(3)])
    # backward decaying solution from a fixed end vector
    wl = _eps(a + 1, io)*al/2; rl = np.sqrt(1 + wl*wl)
    lp, lm = wl + rl, wl - rl
    e_ = np.array([lm if abs(lm) < abs(lp) else lp, 1], dtype=complex); e_ = e_/np.linalg.norm(e_)
    Psi = [None]*(K + 1); dPsi = [None]*(K + 1)
    Psi[K] = e_; dPsi[K] = [np.zeros(2, complex) for _ in range(3)]
    for i in range(K - 1, -1, -1):
        Mi = np.linalg.inv(Ms[i])
        Psi[i] = Mi @ Psi[i + 1]
        dPsi[i] = [Mi @ (dPsi[i + 1][l] - dMs[i][l] @ Psi[i]) for l in range(3)]
    # normalise (constants, chosen at p0)
    F0 = []; F1 = []
    cg = 0 if np.linalg.norm(Phi[K][:, 0]) >= np.linalg.norm(Phi[K][:, 1]) else 1
    for i in range(K + 1):
        if i < i_s:
            mu = np.max(np.linalg.norm(Phi[i], axis=0))
            F0.append(Phi[i]/mu); F1.append([dPhi[i][l]/mu for l in range(3)])
        else:
            ng = np.linalg.norm(Phi[i][:, cg]); nd = np.linalg.norm(Psi[i])
            F0.append(np.column_stack([Phi[i][:, cg]/ng, Psi[i]/nd]))
            F1.append([np.column_stack([dPhi[i][l][:, cg]/ng, dPsi[i][l]/nd]) for l in range(3)])
    return F0, F1, P0


def _am(A):
    return acb_mat([[acb(complex(A[i, j])) for j in range(2)] for i in range(2)])


def _box(F0, F1, D):
    """acb ball matrix F0 + sum_l F1[l] [-D_l, D_l]."""
    out = _am(F0)
    for l in range(3):
        out = out + _am(F1[l])*acb(arb(0, D[l]))
    return out


def lemmaF2_log(s_lo, s_hi, t_lo, t_hi, a_lo, a_hi, i_lo, i_hi, K=200, i_s_list=None, verbose=False, prec=96):
    old = ctx.prec; ctx.prec = prec
    tm, am, im = (t_lo + t_hi)/2, (a_lo + a_hi)/2, (i_lo + i_hi)/2
    D = [(t_hi - t_lo)/2, (a_hi - a_lo)/2, (i_hi - i_lo)/2]
    sm = (s_lo + s_hi)/2
    sigB = arb_iv(s_lo, s_hi); tauB = arb_iv(t_lo, t_hi); aB = arb_iv(a_lo, a_hi); ioB = arb_iv(i_lo, i_hi)
    alB = acb(sigB, tauB)
    al0 = acb(arb(sm), arb(tm))
    I = acb(0, 1)
    # ball and centre matrices, derivative matrices (balls)
    MB = []; M0 = []; dMB = []; lr = []
    for i in range(K):
        mB = aB + (K - i); m0 = arb(am) + (K - i)
        A_ = 2 - (2*mB - 1)*ioB; B_ = 2 - (mB - 1)*ioB
        eB = 2*A_/(mB*B_)
        e0 = 2*(2 - (2*m0 - 1)*arb(im))/(m0*(2 - (m0 - 1)*arb(im)))
        MB.append(acb_mat([[eB*alB, acb(1)], [acb(1), acb(0)]]))
        # centre matrix uses sigma ball (zeroth order in sigma)
        M0.append(acb_mat([[e0*acb(sigB, arb(tm)), acb(1)], [acb(1), acb(0)]]))
        dedio = -4/(B_*B_)
        dedm = 2*(-2*ioB*mB*B_ - A_*(B_ - mB*ioB))/(mB*B_)**2
        z = acb(0)
        dMB.append([acb_mat([[eB*I, z], [z, z]]), acb_mat([[dedm*alB, z], [z, z]]), acb_mat([[dedio*alB, z], [z, z]])])
        rat = 1 - 2*(1 - (mB - 1)*ioB)/(mB*(2 - (mB - 1)*ioB))
        lr.append(float((arb(s_lo)*rat.log()).upper()))
    mA = aB + K
    eA = 2*(2 - (2*mA - 1)*ioB)/(mA*(2 - (mA - 1)*ioB))
    w0 = eA*alB/2; r0 = (1 + w0*w0).sqrt()
    PA = acb_mat([[w0 + r0, w0 - r0], [acb(1), acb(1)]])
    lphi = phase1_log(s_hi, t_lo, t_hi, a_lo, K)
    if i_s_list is None:
        i_s_list = [K - j for j in (10, 20, 30, 45, 60, 80, 110) if j < K]
    best = None
    for i_s in i_s_list:
        F0, F1, _ = build_frames(sm, tm, am, im, K, i_s)
        FB = [_box(F0[i], F1[i], D) for i in range(K + 1)]
        try:
            QB = [f.inv() for f in FB]
        except ZeroDivisionError:
            if verbose: print(f'   i_s={i_s}: singular frame ball', flush=True)
            continue
        Fp = [_am(f) for f in F0]
        Qp = [f.inv() for f in Fp]
        tot = lphi + math.log(_inf_norm(QB[0]*PA))
        worst = 0.0
        for i in range(K):
            E0 = Qp[i + 1]*M0[i]*Fp[i]
            Eb = E0
            for l in range(3):
                Fl1 = _am(F1[i + 1][l]); Fl0 = _am(F1[i][l])
                dE = -(QB[i + 1]*Fl1*QB[i + 1]*MB[i]*FB[i]) + QB[i + 1]*dMB[i][l]*FB[i] + QB[i + 1]*MB[i]*Fl0
                Eb = Eb + dE*acb(arb(0, D[l]))
            nE = _inf_norm(Eb); n0 = _inf_norm(E0)
            worst = max(worst, nE/n0)
            tot += math.log(nE) + lr[i]
        tot += math.log(up(abs(FB[K][0, 0]) + abs(FB[K][0, 1])))
        if verbose:
            print(f'   i_s={i_s}: log bound {tot:.3f} (worst step ratio {worst:.4f})', flush=True)
        if best is None or tot < best[0]:
            best = (tot, i_s)
    ctx.prec = old
    return best


if __name__ == '__main__':
    import time
    L = 4.5; d = math.exp(-L); tq = lambda dd: 1 - math.sqrt(1 - dd)
    io = tq(d)/998.0
    t0 = time.time()
    for K in (200, 400):
        for box in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (4, 1, 1)]:
            da, dt, di = box
            r = lemmaF2_log(2.5, 2.5, 1031.5 - dt/2, 1031.5 + dt/2, 998.0, 998.0 + da, io*(1 - 0.25*di), io*(1 + 0.25*di), K)
            print(K, box, r, f'{time.time()-t0:.1f}s', flush=True)
