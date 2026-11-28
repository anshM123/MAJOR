"""Combined dual SOCP: diagonal rules (endpoint, R_k, M_n) + helicity-flip rules (F_0, F_1, F_3, F_{+-n}), parity-invariant
theory (c = t4/4 real).  Per state the quadratic form in (alpha, beta) (b0 and b4 amplitudes) must be PSD on heavy states and
>= -Theta*(alpha^2 W0 + beta^2 W4) on light states (2<twist<g).  For fixed c0:  maximise Phi = L2 c0^2 - L1 c0  with
3 nu + sum mu = 3; c0 is excluded (with D=0, E'_light=0) iff Phi > 3.  Bisection gives c*, tau* = 4 c*.
State types: (i) parity-even even J>=4: 2x2 PSD with cross term; (ii) parity-odd even J>=4: b0 only, eta=-1;
(iii) J=0,2: b0 only, eta=+-1; (iv) odd J>=5: b4 only."""
import sys, time, numpy as np, cvxpy as cp
from vec import rmat, hmat, Wv
from flip import flip_weights
from core import GAMMA


def head_weights(J, D):
    C = J * (J + 1.0); T = (D - 1) * (D - 2)
    f0 = 1 + (T - 3 * C) / 3.0
    f1 = 2 - 4 * C / 3.0
    return f0, f1


def build(g, K, N, Nf, twist_list, Jlin, Jmax, nJlog):
    Js = np.array(sorted(set(list(range(0, Jlin + 1)) + [int(x) for x in np.geomspace(Jlin + 1, Jmax, nJlog)])), float)
    JJ, TW = np.meshgrid(Js, twist_list, indexing='ij'); J = JJ.ravel(); D = (JJ + TW).ravel()
    out = {}
    for b in (0, 4):
        bb = np.full_like(J, b)
        KK = max(K, N + 1)
        r = rmat(bb, J, D, KK); h = hmat(bb, J, D, N, r) if N >= 2 else np.zeros((len(J), 1))
        w = Wv(bb, D)
        sig = 1.0 if b == 4 else -1.0
        cols = [w] + [h[:, n] for n in range(2, N + 1)] + [sig * r[:, k] for k in range(K + 1)]
        out[b] = (np.array(cols).T, w)
    f0, f1 = head_weights(J, D)
    Fw, rhs, lab = flip_weights(J, D, Nf) if Nf >= 4 else (np.zeros((len(J), 0)), np.zeros(0), [])
    return J, D, out, f0, f1, Fw, rhs


def solve_c0(g, c0, K, N, Nf, theta, Jlin=None, Jmax=None, nJlog=60, nt=24, nl=24, use_flip=True, verbose=False):
    Jlin = Jlin or 6 * g; Jmax = Jmax or 3000 * g
    th_h = g + np.concatenate([[0.0], np.geomspace(1e-3, 1e4, nt)]) * g
    th_l = np.concatenate([[2.0 + 1e-9, 2.05, 2.5, 3, 4, 6], np.geomspace(8, g - 1e-6, nl)])
    nM = max(N - 1, 0); nv = 1 + nM + K + 1
    xs = cp.Variable(nv)              # scaled: x = xs / scx
    J_, D_, out_, f0_, f1_, Fw_, rhs_ = build(g, K, N, Nf, th_h, Jlin, Jmax, nJlog)
    scx = np.maximum(np.abs(out_[0][0] / out_[0][1][:, None]).max(axis=0), np.abs(out_[4][0] / out_[4][1][:, None]).max(axis=0)); scx[scx == 0] = 1
    x = cp.multiply(1.0 / scx, xs)
    nf = 2 + (2 * (Nf - 3) - 1 + 1 if Nf >= 4 else 0)
    nk = 2 + (1 + 2 * (Nf - 3) if Nf >= 4 else 0)
    sw = np.sqrt(Wv(np.zeros_like(J_), D_) * Wv(4 * np.ones_like(J_), D_))
    sck = np.ones(nk)
    sck[0] = np.abs(f0_ / Wv(np.zeros_like(J_), D_)).max(); sck[1] = np.abs(f1_ / Wv(np.zeros_like(J_), D_)).max()
    if Nf >= 4:
        sck[2:] = np.abs(Fw_ / sw[:, None]).max(axis=0)
    sck[sck == 0] = 1
    ks = cp.Variable(nk)
    kap = cp.multiply(1.0 / sck, ks)
    cons = [xs[0] >= 0, 3 * x[0] + cp.sum(x[1:1 + nM]) == 3]
    for twl, light in ((th_h, False), (th_l, True)):
        J, D, out, f0, f1, Fw, rhs = build(g, K, N, Nf, twl, Jlin, Jmax, nJlog)
        A0, W0 = out[0]; A4, W4 = out[4]
        s0 = W0; s4 = W4; sx = np.sqrt(W0 * W4)
        M0 = (A0 / s0[:, None]) / scx[None, :]; M4 = (A4 / s4[:, None]) / scx[None, :]
        H0 = np.zeros((len(J), nk)); H0[:, 0] = f0 / s0 / sck[0]; H0[:, 1] = f1 / s0 / sck[1]
        F0p = M0 @ xs + H0 @ ks
        F0m = M0 @ xs - H0 @ ks
        F4 = M4 @ xs
        sl = theta if light else 0.0
        even = (J % 2 == 0)
        big = even & (J >= 4); small = even & (J < 4); odd = (J % 2 == 1) & (J >= 5)
        # (iii) J=0,2 both parities ; (ii) parity-odd even J>=4 ; (iv) odd J
        cons += [F0p[np.where(small)[0]] >= -sl, F0m[np.where(small | big)[0]] >= -sl, F4[np.where(odd)[0]] >= -sl]
        ib = np.where(big)[0]
        if use_flip and Nf >= 4:
            G = (Fw[ib] / sx[ib, None] / sck[None, 2:]) @ ks[2:]
        else:
            G = np.zeros(len(ib))
        a = F0p[ib] + sl; bq = F4[ib] + sl
        cons += [cp.SOC(a + bq, cp.vstack([a - bq, G]))]
    L2 = sum(GAMMA.get(n, 0.0) * x[1 + n - 2] for n in range(2, N + 1)) - x[1 + nM + 2]
    L1 = kap[0] + 2 * kap[1]
    if Nf >= 4:
        L1 = L1 - 2 * kap[2] + cp.sum(kap[3:])
    if not use_flip:
        cons += [ks == 0]
    prob = cp.Problem(cp.Maximize(L2 * c0 ** 2 - L1 * c0), cons)
    try:
        prob.solve(solver='CLARABEL')
    except Exception as e:
        return None, 'error'
    return prob.value, prob.status


def bisect(g, K, N, Nf, theta, use_flip=True, lo=1e-6, hi=1.0, it=22, **kw):
    for i in range(it):
        mid = np.sqrt(lo * hi)
        v, st = solve_c0(g, mid, K, N, Nf, theta, use_flip=use_flip, **kw)
        if v is None or st not in ('optimal', 'optimal_inaccurate'):
            v = np.inf if st == 'unbounded' else -1
        if v > 3: hi = mid
        else: lo = mid
    return hi


if __name__ == '__main__':
    g = int(sys.argv[1]); p = float(sys.argv[2]); K = int(sys.argv[3]); N = int(sys.argv[4]); Nf = int(sys.argv[5])
    th = g ** p
    for uf in (False, True):
        t0 = time.time()
        cs = bisect(g, K, N, Nf, th, use_flip=uf)
        print(f"g={g} Theta=g^{p} K={K} N={N} Nf={Nf} flip={uf}: c*={cs:.4e} tau*={4*cs:.4e} tau*g={4*cs*g:.3f} [{time.time()-t0:.0f}s]", flush=True)
