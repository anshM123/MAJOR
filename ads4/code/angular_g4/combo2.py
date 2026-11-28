"""Combined dual (diagonal + flip) as an LP via a per-state inner approximation of the 2x2 PSD condition:
|G^| <= 2 F0^/s  and  |G^| <= 2 s F4^   (sufficient for G^2 <= 4 F0 F4; s updated iteratively, each LP is a valid bound).
Hats: F0^ = F0/W0, F4^ = F4/W4, G^ = G/sqrt(W0 W4).  Light states: F^ >= -Theta (form + Theta*(a^2 W0 + b^2 W4) PSD).
For fixed c0: maximise Phi = L2 c0^2 - L1 c0 (3 nu + sum mu = 3); c0 excluded iff Phi > 3 (D = 0, E'_light = 0)."""
import sys, time, numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack as svstack, hstack as shstack
from vec import rmat, hmat, Wv
from flip import flip_weights
from core import GAMMA
from combo import head_weights


def blocks(g, K, N, Nf, tw, Jlin, Jmax, nJlog):
    Js = np.array(sorted(set(list(range(0, Jlin + 1)) + [int(x) for x in np.geomspace(Jlin + 1, Jmax, nJlog)])), float)
    JJ, TW = np.meshgrid(Js, tw, indexing='ij'); J = JJ.ravel(); D = (JJ + TW).ravel()
    M = {}
    for b in (0, 4):
        bb = np.full_like(J, b); KK = max(K, N + 1)
        r = rmat(bb, J, D, KK); h = hmat(bb, J, D, N, r) if N >= 2 else np.zeros((len(J), 1))
        w = Wv(bb, D); sig = 1.0 if b == 4 else -1.0
        M[b] = np.array([w] + [h[:, n] for n in range(2, N + 1)] + [sig * r[:, k] for k in range(K + 1)]).T / w[:, None]
    f0, f1 = head_weights(J, D); W0 = Wv(np.zeros_like(J), D)
    H = np.stack([f0 / W0, f1 / W0], axis=1)
    if Nf >= 4:
        Fw, rhs, lab = flip_weights(J, D, Nf); G = Fw / np.sqrt(W0 * Wv(4 * np.ones_like(J), D))[:, None]
    else:
        G = np.zeros((len(J), 0))
    return J, D, M, H, G


class Combo:
    def __init__(self, g, K, N, Nf, theta, use_flip=True, Jlin=None, Jmax=None, nJlog=60, nt=30, nl=30):
        self.g, self.K, self.N, self.Nf, self.theta, self.uf = g, K, N, Nf, theta, use_flip
        Jlin = Jlin or 6 * g; Jmax = Jmax or 3000 * g
        th_h = g + np.concatenate([[0.0], np.geomspace(1e-3, 1e4, nt)]) * g
        self.parts = [(blocks(g, K, N, Nf, th_h, Jlin, Jmax, nJlog), 0.0)]
        if theta is not None:
            th_l = np.concatenate([[2.0 + 1e-9, 2.05, 2.5, 3, 4, 6], np.geomspace(8, g - 1e-6, nl)])
            self.parts.append((blocks(g, K, N, Nf, th_l, Jlin, Jmax, nJlog), theta))
        self.nx = 1 + max(N - 1, 0) + K + 1
        self.nk = 2 + ((1 + 2 * (Nf - 3)) if (Nf >= 4 and use_flip) else 0)
        if not use_flip: self.nk = 0
        self.s = [None for _ in self.parts]

    def lp(self, c0):
        nx, nk = self.nx, self.nk; nM = max(self.N - 1, 0)
        nv = nx + nk
        rows, rhs = [], []
        for ip, ((J, D, M, H, G), sl) in enumerate(self.parts):
            ev = (J % 2 == 0); big = ev & (J >= 4); small = ev & (J < 4); odd = (J % 2 == 1) & (J >= 5)
            Hk = np.zeros((len(J), nk))
            if nk: Hk[:, :2] = H
            Gk = np.zeros((len(J), nk))
            if nk > 2: Gk[:, 2:] = G
            A0p = np.hstack([M[0], Hk]); A0m = np.hstack([M[0], -Hk]); A4 = np.hstack([M[4], np.zeros((len(J), nk))])
            # F >= -sl   ->  -F <= sl
            for A, mask in ((A0p, small), (A0m, small | big), (A4, odd)):
                rows.append(-A[mask]); rhs.append(np.full(mask.sum(), sl))
            ib = np.where(big)[0]
            s = self.s[ip][ib] if self.s[ip] is not None else np.ones(len(ib))
            Gb = np.hstack([np.zeros((len(ib), nx)), Gk[ib][:, :nk]]) if nk else np.zeros((len(ib), nv))
            F0b = A0p[ib]; F4b = A4[ib]
            # |G| <= 2 (F0 + sl)/s  and |G| <= 2 s (F4 + sl)
            for sg in (1, -1):
                rows.append(sg * Gb - 2 * F0b / s[:, None]); rhs.append(2 * sl / s)
                rows.append(sg * Gb - 2 * s[:, None] * F4b); rhs.append(2 * sl * s)
        Aub = np.vstack(rows); bub = np.concatenate(rhs)
        sc = np.abs(Aub).max(axis=0); sc[sc == 0] = 1
        Aeq = np.zeros((1, nv)); Aeq[0, 0] = 3; Aeq[0, 1:1 + nM] = 1
        cobj = np.zeros(nv)
        for n in range(2, self.N + 1): cobj[1 + n - 2] = GAMMA.get(n, 0.0) * c0 ** 2
        cobj[1 + nM + 2] = -c0 ** 2
        if nk:
            l1 = np.zeros(nk); l1[0] = 1; l1[1] = 2
            if nk > 2: l1[2] = -2; l1[3:] = 1
            cobj[nx:] = -l1 * c0
        bounds = [(0, None)] + [(None, None)] * (nv - 1)
        res = linprog(-cobj / sc, A_ub=Aub / sc, b_ub=bub, A_eq=Aeq / sc, b_eq=[3.0], bounds=bounds, method='highs')
        if res.status == 3: return np.inf, None
        if res.status != 0: return None, None
        return -res.fun, res.x / sc

    def update_s(self, x):
        nx, nk = self.nx, self.nk
        for ip, ((J, D, M, H, G), sl) in enumerate(self.parts):
            Hk = np.zeros((len(J), nk));
            if nk: Hk[:, :2] = H
            F0 = np.hstack([M[0], Hk]) @ x + sl; F4 = np.hstack([M[4], np.zeros((len(J), nk))]) @ x + sl
            self.s[ip] = np.sqrt(np.clip(F0, 1e-300, None) / np.clip(F4, 1e-300, None))

    def cstar(self, lo, hi, it=18, iters_s=3):
        for i in range(it):
            mid = np.sqrt(lo * hi)
            self.s = [None for _ in self.parts]
            v, x = self.lp(mid)
            for _ in range(iters_s):
                if v is None or not np.isfinite(v) or v > 3: break
                self.update_s(x); v, x = self.lp(mid)
            if v is not None and v > 3: hi = mid
            else: lo = mid
        return hi


if __name__ == '__main__':
    g = int(sys.argv[1]); p = None if sys.argv[2] == 'inf' else float(sys.argv[2])
    K, N, Nf = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    th = None if p is None else g ** p
    for uf in (False, True):
        t0 = time.time()
        cm = Combo(g, K, N, Nf, th, use_flip=uf)
        cs = cm.cstar(1e-5, 1.0)
        print(f"g={g} p={p} K={K} N={N} Nf={Nf} flip={uf}: tau*={4*cs:.4e} tau*g={4*cs*g:.3f} [{time.time()-t0:.0f}s]", flush=True)
