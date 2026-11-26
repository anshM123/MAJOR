"""
tpf3.py: Lemma F (turning point) with first-order Riccati frames and a second-order defect bound.  ANALYTIC_SLIVER2 §3.

Box: p = (tau, a, iota, sigma) in B = prod [p0_l - D_l, p0_l + D_l] (sigma >= 0, iota = 1/N), b = 0 (s = 1).
Last segment i = 0..K-1 <-> k = k_A + i, m = N - k = a + K - i;  eps = 2/m - 2 iota/(2 - (m-1) iota) (= eps_k exactly),
w_i = eps alpha/2, M_i = [[2 w_i, 1], [1, 0]].
Frames: F_i(p) = [[R1_i(p), R2_i(p)], [1, 1]] with AFFINE ratios R_c(p) = R_c^0 + sum_l g_{c,l} (p_l - p0_l) (point coefficients from the
Riccati recursion R_{i+1} = 2 w_i + 1/R_i and its p-derivative at p0; segment 1: two forward ratios from lam+-(M_{k_A}(p0));
segment 2 (i >= i_s): forward ratio 1 and the backward ratio R_i = 1/(R_{i+1} - 2 w_i) from a fixed end value).
Exact identity (for i + 1 != i_s):  M_i F_i = F_{i+1} diag(R1_i, R2_i) + [[R1_i delta1_i, R2_i delta2_i], [0, 0]],
    delta_c(p) := 2 w_i(p) + 1/R_{c,i}(p) - R_{c,i+1}(p)      (the Riccati defect of the affine frames).
Hence E_i = Q_{i+1} M_i F_i = diag(R1_i, R2_i) + Q_{i+1} [[R1 delta1, R2 delta2], [0, 0]] and
    ||E_i||_inf <= max_c |R_c(p)| + max_r |(Q_{i+1})_{r,1}| (|R1 delta1| + |R2 delta2|),
    |R(p)| <= sqrt(|R^0|^2 + 2 sum_l |Re(conj(R^0) g_l)| D_l + (sum_l |g_l| D_l)^2),
    |delta(p)| <= |delta(p0)| + sum_l |d_l delta(p0)| D_l + (1/2) sum_{l,m} sup_B |d_l d_m delta| D_l D_m,
    d_l delta = 2 d_l w - g_l/R^2 - g'_l,   d_l d_m delta = 2 d_l d_m w + 2 g_l g_m / R^3   (R over B as an Arb ball).
The switch step and the boundary factors are evaluated as plain Arb balls over B.  Phase I: tpframes.phase1_log.
    |S^_n(alpha)| <= Phi_I ||Q_0 P_{k_A}(B)||_inf prod_i [(d_{k+1}/d_k)^{sigma_lo} ||E_i||_inf] (|R1_K| + |R2_K|).
"""
import math, sys
import numpy as np
sys.path.insert(0, '.')
from flint import arb, acb, ctx, acb_mat
from tpframes import phase1_log, up, arb_iv, _inf_norm
ctx.prec = 96
NP = 4   # parameters: 0 tau, 1 a, 2 iota, 3 sigma


def eps_derivs(m, io):
    """eps and its derivatives in (m, iota): returns dict (works for floats or arb)."""
    B = 2 - (m - 1)*io
    e = 2/m - 2*io/B
    return dict(e=e, m=-2/m**2 - 2*io*io/B**2, i=-4/B**2, mm=4/m**3 - 4*io**3/B**3, mi=-8*io/B**3, ii=-8*(m - 1)/B**3)


def w_derivs(ed, al):
    """first and second derivatives of w = eps alpha/2 in (tau, a, iota, sigma); alpha = sigma + i tau, dm/da = 1."""
    I = 1j if not isinstance(al, acb) else acb(0, 1)
    e = ed['e']
    d1 = [I*e/2, ed['m']*al/2, ed['i']*al/2, e/2]
    z = 0*e
    d2 = [[z, I*ed['m']/2, I*ed['i']/2, z],
          [I*ed['m']/2, ed['mm']*al/2, ed['mi']*al/2, ed['m']/2],
          [I*ed['i']/2, ed['mi']*al/2, ed['ii']*al/2, ed['i']/2],
          [z, ed['m']/2, ed['i']/2, z]]
    return d1, d2


def build(sig, tau, a, io, K, i_s):
    al = complex(sig, tau)
    w = []; dw = []
    for i in range(K):
        ed = eps_derivs(a + K - i, io)
        w.append(ed['e']*al/2)
        dw.append(w_derivs(ed, al)[0])
    w0 = eps_derivs(a + K, io)['e']*al/2; r0 = np.sqrt(1 + w0*w0)
    R = [[w0 + r0, w0 - r0]]; g = [[[0j]*NP, [0j]*NP]]
    for i in range(K):
        Rn = []; gn = []
        for c in range(2):
            Rc = R[i][c]
            Rn.append(2*w[i] + 1/Rc)
            gn.append([2*dw[i][l] - g[i][c][l]/Rc**2 for l in range(NP)])
        R.append(Rn); g.append(gn)
    wl = w[K - 1]; rl = np.sqrt(1 + wl*wl)
    lp, lm = wl + rl, wl - rl
    Rd = [None]*(K + 1); gd = [None]*(K + 1)
    Rd[K] = lm if abs(lm) < abs(lp) else lp; gd[K] = [0j]*NP
    for i in range(K - 1, -1, -1):
        Rd[i] = 1/(Rd[i + 1] - 2*w[i])
        gd[i] = [-Rd[i]**2*(gd[i + 1][l] - 2*dw[i][l]) for l in range(NP)]
    cols = []
    for i in range(K + 1):
        if i < i_s:
            cols.append(((R[i][0], g[i][0]), (R[i][1], g[i][1])))
        else:
            cols.append(((R[i][0], g[i][0]), (Rd[i], gd[i])))
    return cols


def _c(z):
    return acb(complex(z))


def lemmaF3_log(s_lo, s_hi, t_lo, t_hi, a_lo, a_hi, i_lo, i_hi, K=200, i_s_list=None, verbose=False, prec=96):
    old = ctx.prec; ctx.prec = prec
    p0 = [(t_lo + t_hi)/2, (a_lo + a_hi)/2, (i_lo + i_hi)/2, (s_lo + s_hi)/2]
    D = [(t_hi - t_lo)/2, (a_hi - a_lo)/2, (i_hi - i_lo)/2, (s_hi - s_lo)/2]
    Da = [arb(x) for x in D]
    box = [arb_iv(t_lo, t_hi), arb_iv(a_lo, a_hi), arb_iv(i_lo, i_hi), arb_iv(s_lo, s_hi)]
    alB = acb(box[3], box[0]); al0 = acb(arb(p0[3]), arb(p0[0]))
    # per-step Arb data: w at p0 (point), dw at p0, second derivatives over B, w over B, d-ratio
    W0 = []; DW0 = []; D2B = []; WB = []; lr = []
    for i in range(K):
        m0 = arb(p0[1]) + (K - i); mB = box[1] + (K - i)
        ed0 = eps_derivs(m0, arb(p0[2])); edB = eps_derivs(mB, box[2])
        d10, _ = w_derivs(ed0, al0)
        _, d2B = w_derivs(edB, alB)
        W0.append(ed0['e']*al0/2); DW0.append(d10); D2B.append(d2B)
        WB.append(edB['e']*alB/2)
        rat = 1 - 2*(1 - (mB - 1)*box[2])/(mB*(2 - (mB - 1)*box[2]))
        lr.append(float((arb(s_lo)*rat.log()).upper()))
    mA = box[1] + K
    wA = eps_derivs(mA, box[2])['e']*alB/2; rA = (1 + wA*wA).sqrt()
    PA = acb_mat([[wA + rA, wA - rA], [acb(1), acb(1)]])
    lphi = phase1_log(s_hi, t_lo, t_hi, a_lo, K)
    if i_s_list is None:
        i_s_list = [K - j for j in (10, 20, 30, 45, 60, 80) if j < K]
    best = None
    for i_s in i_s_list:
        cols = build(p0[3], p0[0], p0[1], p0[2], K, i_s)
        # affine frame balls and the |R| bounds
        def Rball(cg):
            R0, gg = cg
            v = _c(R0)
            for l in range(NP):
                v = v + _c(gg[l])*acb(arb(0, D[l]))
            return v
        def Rabs(cg):
            R0, gg = cg
            R0a = _c(R0)
            lin = arb(0); tot1 = arb(0)
            for l in range(NP):
                lin += abs((R0a.conjugate()*_c(gg[l])).real)*Da[l]
                tot1 += abs(_c(gg[l]))*Da[l]
            return (abs(R0a)**2 + 2*lin + tot1**2).sqrt()
        FB = [acb_mat([[Rball(cols[i][0]), Rball(cols[i][1])], [acb(1), acb(1)]]) for i in range(K + 1)]
        try:
            QB = [f.inv() for f in FB]
        except ZeroDivisionError:
            if verbose: print(f'   i_s={i_s}: singular', flush=True)
            continue
        F0p = acb_mat([[_c(cols[0][0][0]), _c(cols[0][1][0])], [acb(1), acb(1)]])
        tot = lphi + math.log(_inf_norm(F0p.inv()*PA))
        for i in range(K):
            if i + 1 == i_s:
                E = QB[i + 1]*acb_mat([[2*WB[i], acb(1)], [acb(1), acb(0)]])*FB[i]
                tot += math.log(_inf_norm(E)) + lr[i]
                continue
            dmax = arb(0); rsum = arb(0)
            Rmax = arb(0)
            for c in range(2):
                (R0, gl), (R1, gl1) = cols[i][c], cols[i + 1][c]
                R0a, R1a = _c(R0), _c(R1)
                d0 = 2*W0[i] + 1/R0a - R1a
                bound = abs(d0)
                for l in range(NP):
                    dl = 2*DW0[i][l] - _c(gl[l])/R0a**2 - _c(gl1[l])
                    bound += abs(dl)*Da[l]
                RB = Rball(cols[i][c])
                for l in range(NP):
                    for mm in range(NP):
                        if D[l] == 0 or D[mm] == 0:
                            continue
                        h = 2*D2B[i][l][mm] + 2*_c(gl[l])*_c(gl[mm])/RB**3
                        bound += abs(h)*Da[l]*Da[mm]/2
                Ra = Rabs(cols[i][c])
                Rmax = Rmax if Rmax > Ra else arb(up(Ra)) if Rmax < Ra else arb(max(up(Rmax), up(Ra)))
                rsum += Ra*bound
            qrow = max(up(abs(QB[i + 1][0, 0])), up(abs(QB[i + 1][1, 0])))
            nE = up(Rmax) + qrow*up(rsum)
            tot += math.log(nE) + lr[i]
        tot += math.log(up(Rabs(cols[K][0]) + Rabs(cols[K][1])))
        if verbose:
            print(f'   i_s={i_s}: log bound {tot:.3f}', flush=True)
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
        for box in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (4, 1, 1), (10, 1, 1)]:
            da, dt, di = box
            r = lemmaF3_log(2.5, 2.5, 1031.5 - dt/2, 1031.5 + dt/2, 998.0, 998.0 + da, io*(1 - 0.25*di), io*(1 + 0.25*di), K)
            print(K, box, r, f'{time.time()-t0:.1f}s', flush=True)
