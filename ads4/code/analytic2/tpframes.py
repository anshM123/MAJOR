"""
tpframes.py: Lemma F (turning point by exact-solution frames) -- certified bound of |S^_n(alpha)|, alpha = sigma + i tau, b = 0 (s = 1),
uniformly over a box  sigma in [s_lo, s_hi] (s_lo >= 0), tau in [t_lo, t_hi], a in [a_lo, a_hi], iota = 1/N in [i_lo, i_hi],
every integer n > K.  ANALYTIC_SLIVER2 §3.

Recurrence (dual Hahn, s = 1): v_{k+1} = M_k v_k, v_k = (S_k, S_{k-1}), M_k = [[2 w_k, 1], [1, 0]], w_k = eps_k alpha/2,
eps_k = 2(2k+1)/((N-k)(N+k+1)).  In the variable m = N - k:
    eps = 2(2 - (2m-1) iota)/(m (2 - (m-1) iota)),     d_{k+1}/d_k = 1 - 2(1 - (m-1) iota)/(m (2 - (m-1) iota))   (d_k = 1 - k(k+1)/(N(N+1))).
Phase I (k < k_A := n - K, i.e. m > a + K): Theorem E' (ANALYTIC_SLIVER §5.2) stopped at k_A, with a' = a + K:
    d_{k_A}^sigma |c_{k_A}|_inf <= Phi_I := c0 (2/(1+sqrt(1-k'^2)))^sigma exp(q[asin k' + (q/2) log(1/(1-k'^2))]),  k' = tau/(a+K),
    c_k = P_k^{-1} v_k eigen-coordinates, P_k = [[lam+, lam-], [1, 1]], q = |alpha|/tau, c0 = max(|1-lam-|, |lam+-1|)/(2|r|) at k = 0.
Last segment (k = k_A .. n-1, m = a+K .. a+1): any point frames F_k (k = k_A..n), Q_k = F_k^{-1} (Arb):
    |S^_n| <= Phi_I ||Q_{k_A} P_{k_A}||_inf  prod_k [ (d_{k+1}/d_k)^sigma ||Q_{k+1} M_k F_k||_inf ]  ||row_1 F_n||_1,
M_k, P_{k_A} as Arb balls over the box.  Frames: two exact solutions at the box centre (traveling waves started from the
eigenvectors of M_{k_A}, common normalisation), switched at index i_s to (forward solution, backward-decaying solution) with
separate normalisations.  Rigorous for any frames; the choice only affects sharpness.
"""
import math, sys
import numpy as np
from flint import arb, acb, ctx, acb_mat
ctx.prec = 96


def up(x):
    return float(arb(x).upper())


def arb_iv(lo, hi):
    lo = arb(lo); hi = arb(hi)
    return arb((lo + hi)/2, ((hi - lo)/2).upper()*(1 + 1e-12) + 1e-300)


def _inf_norm(Mm):
    return max(up(abs(Mm[i, 0]) + abs(Mm[i, 1])) for i in range(2))


def phase1_log(s_hi, t_lo, t_hi, a_lo, K):
    """log Phi_I (Arb upper bound), monotone: kappa' = t_hi/(a_lo+K), q = |alpha|/tau <= sqrt(s_hi^2 + t_lo^2)/t_lo ... use max."""
    s = arb(s_hi); th = arb(t_hi); tl = arb(t_lo)
    kp = th/(arb(a_lo) + K)
    assert kp < 1
    q = (s*s + tl*tl).sqrt()/tl          # |alpha|/tau decreasing in tau
    sq = (1 - kp*kp).sqrt()
    al = (s*s + th*th).sqrt()
    c0 = 1 + 2*al/arb(a_lo)**2           # c0 = 1 + O(|alpha|/N^2) <= 1 + 2|alpha|/a^2 (E' statement)
    v = c0.log() + s*(2/(1 + sq)).log() + q*(kp.asin() + q/2*(-(1 - kp*kp).log()))
    return up(v)


def _eps_np(m, io):
    return 2*(2 - (2*m - 1)*io)/(m*(2 - (m - 1)*io))


def frames(sig, tau, a, io, K, i_s):
    """point frames F_k, k = 0..K (index i <-> k = k_A + i, m = a + K - i), at a single parameter point."""
    al = complex(sig, tau)
    Ms = []
    for i in range(K):
        m = a + K - i
        w = _eps_np(m, io)*al/2
        Ms.append(np.array([[2*w, 1], [1, 0]], dtype=complex))
    w0 = _eps_np(a + K, io)*al/2
    r0 = np.sqrt(1 + w0*w0)
    P0 = np.array([[w0 + r0, w0 - r0], [1, 1]], dtype=complex)
    # forward pair
    Phi = [P0.copy()]
    for i in range(K):
        Phi.append(Ms[i] @ Phi[-1])
        nr = np.max(np.linalg.norm(Phi[-1], axis=0))
        Phi[-1] = Phi[-1]/nr
    # backward decaying solution from the end
    wl = _eps_np(a + 1, io)*al/2
    rl = np.sqrt(1 + wl*wl)
    lp, lm = wl + rl, wl - rl
    vend = np.array([lm if abs(lm) < abs(lp) else lp, 1], dtype=complex)
    back = [None]*(K + 1); back[K] = vend/np.linalg.norm(vend)
    for i in range(K - 1, -1, -1):
        v = np.linalg.solve(Ms[i], back[i + 1])
        back[i] = v/np.linalg.norm(v)
    F = []
    for i in range(K + 1):
        if i < i_s:
            mu = np.max(np.linalg.norm(Phi[i], axis=0))
            F.append(Phi[i]/mu)
        else:
            g = Phi[i][:, 0] + 0*Phi[i][:, 1]
            # forward growing direction: the larger of the two forward columns
            c0n = np.linalg.norm(Phi[i][:, 0]); c1n = np.linalg.norm(Phi[i][:, 1])
            g = Phi[i][:, 0] if c0n >= c1n else Phi[i][:, 1]
            F.append(np.column_stack([g/np.linalg.norm(g), back[i]]))
    return F, P0


def _acbm(A):
    return acb_mat([[acb(complex(A[i, j])) for j in range(2)] for i in range(2)])


def lemmaF_log(s_lo, s_hi, t_lo, t_hi, a_lo, a_hi, i_lo, i_hi, K=400, i_s_list=None, verbose=False):
    """log upper bound of |S^_n(sigma + i tau)| over the box (b = 0).  Returns min over switch points."""
    sig = arb_iv(s_lo, s_hi); tau = arb_iv(t_lo, t_hi); a = arb_iv(a_lo, a_hi); io = arb_iv(i_lo, i_hi)
    al = acb(sig, tau)
    sm, tm, am, im = (s_lo + s_hi)/2, (t_lo + t_hi)/2, (a_lo + a_hi)/2, (i_lo + i_hi)/2
    # Arb matrices over the box
    Mb = []; lr = []
    for i in range(K):
        m = a + (K - i)
        eps = 2*(2 - (2*m - 1)*io)/(m*(2 - (m - 1)*io))
        w = eps*al/2
        Mb.append(acb_mat([[2*w, acb(1)], [acb(1), acb(0)]]))
        rat = 1 - 2*(1 - (m - 1)*io)/(m*(2 - (m - 1)*io))
        lr.append(rat.log())            # log(d_{k+1}/d_k) < 0
    m = a + K
    eps = 2*(2 - (2*m - 1)*io)/(m*(2 - (m - 1)*io))
    w0 = eps*al/2; r0 = (1 + w0*w0).sqrt()
    PA = acb_mat([[w0 + r0, w0 - r0], [acb(1), acb(1)]])
    lphi = phase1_log(s_hi, t_lo, t_hi, a_lo, K)
    if i_s_list is None:
        i_s_list = [K - j for j in (5, 10, 20, 30, 45, 60, 80, 110, 150, 200) if j < K]
    best = None
    for i_s in i_s_list:
        F, _ = frames(sm, tm, am, im, K, i_s)
        Fa = [_acbm(f) for f in F]
        Qa = [f.inv() for f in Fa]
        tot = lphi + math.log(_inf_norm(Qa[0]*PA))
        for i in range(K):
            E = Qa[i + 1]*Mb[i]*Fa[i]
            tot += math.log(_inf_norm(E)) + up(sig*lr[i]) if False else math.log(_inf_norm(E)) + float((arb(s_lo)*lr[i]).upper())
        tot += math.log(up(abs(Fa[K][0, 0]) + abs(Fa[K][0, 1])))
        if verbose:
            print(f'   i_s={i_s}: log bound {tot:.3f}', flush=True)
        if best is None or tot < best[0]:
            best = (tot, i_s)
    return best


if __name__ == '__main__':
    import time
    t0 = time.time()
    L = 4.5; d = math.exp(-L)
    # iota = 1/N with N = a/t, t ~ d/2
    a_lo, a_hi = 998.0, 999.0
    tq = lambda dd: (1 - math.sqrt(1 - dd))
    i_lo = tq(math.exp(-5.0))/a_hi; i_hi = tq(d)/a_lo
    for K in (200, 400):
        print('K', K, lemmaF_log(2.5, 2.5, 1031.0, 1032.0, a_lo, a_hi, i_lo, i_hi, K, verbose=True), f'{time.time()-t0:.1f}s')
