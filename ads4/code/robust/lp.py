"""LP for a robust flat functional for |ghat_3|^2 (ROBUST_FUNCTIONAL.md section 2).  Floating point: results are [NE].

Variables v = (x_0..x_d2, y_0..y_d3, a):  psi2 = (1-p)^3 p sum x_i p^i,  psi3 = (1-p)^3 sum y_j p^j,  a >= 0.
Normalisation k = -2 mu2_6 + 2 mu3_4 + 4a = 1, so A = 16 x_0, B = 16 Lfin, B_AdS = B - A log(3/2).

Robustness criterion with parameters (delta, C, u0, dsc):
  (R1) physical states J allowed, u = 1/m^2 in (0, u0]:  FF - (delta + C u) N >= 0,   N = |B2| + |B3| + |FW|;
  (R2) scaling limit u -> 0:  2 int psi2 J_0(bhat p) >= dsc * rho(bhat),  rho = rho0 (1 + (bhat/b0)^2)^{-3/2}
       (rho0 = 2 int psi2^{F_g3}/k^{F_g3} = 0.4489, b0 = (80/(2 int psi2^{F_g3}))^{1/3} = 1.865: the F_g3 scale; tail 80/(k_0 bhat^3));
  (R3) large bhat (Lemma M asymptotics, bhat in [BAS, 1e4]):  the two-term asymptotic action satisfies (R1);
  (M)  matter (FLAT_GRAVITON_BOUND.md section 4.2): coefficients of a0 = |g_0|^2 and a2 = |g_2|^2 in the low side <= 0 for 0 < m_l <= 1.
R1 is linear for each sign pattern (s2, s3) of (B2, B3); handled by cutting planes (add the violated pattern).
Objective: minimise A log(LREF) + B_AdS  (LREF = inf: minimise A, tie-break B).
"""
import sys, os, time, json
import numpy as np
from scipy.optimize import linprog
from rcommon import (basis2, basis3, mom, flat_parts_basis, phi_coef, jmax_for, scaling_B2, Functional, SIGMA, C4,
                     peval, paper_fg3)

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache')

RHO0 = 2 * 6.169047619047619 / 27.48509823509836
B0 = (80 / (2 * 6.169047619047619)) ** (1 / 3)
def rho(bh):
    return RHO0 / (1 + (bh / B0) ** 2) ** 1.5

def u_grid(nu=70):
    u = np.unique(np.concatenate([np.logspace(np.log10(2.5e-5), 0, nu), np.linspace(0.5, 1, 26), np.linspace(0.8, 1, 21)]))
    return u

class Grid:
    """All rows: for each (u, h, J): B2 basis row, B3 basis row, FW coefficient (per unit a), u."""
    def __init__(self, d2, d3, r=3, bmax=200.0, Jcap=2500, nu=70, umax=1.0, verbose=True):
        self.d2, self.d3, self.r = d2, d3, r
        fn = os.path.join(CACHE, 'grid_d%d_%d_r%d_b%d_J%d_n%d.npz' % (d2, d3, r, bmax, Jcap, nu))
        if os.path.exists(fn):
            z = np.load(fn)
            self.B2, self.B3, self.FWc, self.u, self.J, self.h = [z[k] for k in ['B2', 'B3', 'FWc', 'u', 'J', 'h']]
        else:
            b2 = basis2(d2, r); b3 = basis3(d3, r)
            B2l, B3l, FWl, ul, Jl, hl = [], [], [], [], [], []
            t0 = time.time()
            for u in u_grid(nu):
                Jm = jmax_for(u, bmax, Jcap)
                for h in ['++', '+-']:
                    Js, B2, B3 = flat_parts_basis(h, u, Jm, b2, b3)
                    if h == '++':
                        sel = Js % 2 == 0
                        Js, B2, B3 = Js[sel], B2[sel], B3[sel]
                    B2l.append(B2); B3l.append(B3); FWl.append(phi_coef(h, Js, u))
                    ul.append(np.full(len(Js), u)); Jl.append(Js); hl.append(np.full(len(Js), SIGMA[h]))
                if verbose:
                    print('  grid u=%.3e Jmax=%d (%.0fs)' % (u, Jm, time.time() - t0), flush=True)
            self.B2 = np.concatenate(B2l); self.B3 = np.concatenate(B3l); self.FWc = np.concatenate(FWl)
            self.u = np.concatenate(ul); self.J = np.concatenate(Jl); self.h = np.concatenate(hl)
            np.savez(fn, B2=self.B2, B3=self.B3, FWc=self.FWc, u=self.u, J=self.J, h=self.h)
        self.n2 = d2 + 1; self.n3 = d3 + 1; self.nv = self.n2 + self.n3 + 1
        # scaling limit
        self.bsc = np.unique(np.concatenate([np.linspace(0, 60, 1201), np.linspace(60, 400, 1700)]))
        self.Bsc = scaling_B2(basis2(d2, r), self.bsc)
        # asymptotic (Lemma M) rows: bhat in [BAS, 1e4], u in (0, 1/4]
        self.BAS = bmax * 0.75
        bb = np.logspace(np.log10(self.BAS), 4, 120)
        uu = np.logspace(-9, np.log10(0.25), 160)
        Bb, Uu = np.meshgrid(bb, uu)
        self.ab, self.au = Bb.ravel(), Uu.ravel()
        b2 = basis2(d2, r); b3 = basis3(d3, r)
        d2pp = np.array([2 * c[2] for c in b2])                     # psi2''(0) per basis
        p3_0 = np.array([c[0] for c in b3]); p3pp = np.array([2 * c[2] if len(c) > 2 else 0. for c in b3])
        self.asym = {}
        for h, s in SIGMA.items():
            A2 = np.outer(1 / self.ab ** 3, -d2pp)                   # B2 ~ -psi2''(0)/b^3
            A3 = s * self.au[:, None] * (np.outer(1 / self.ab, p3_0) - np.outer(1 / (2 * self.ab ** 3), p3pp))
            Fc = self.ab ** 2 * self.au ** 2 / 2 - self.au ** 3 / 2 - (C4[h] - s) * self.au ** 3
            self.asym[h] = (A2, A3, Fc)
        # low side
        self.mu2_6 = np.array([mom(c, 6) for c in b2]); self.mu3_4 = np.array([mom(c, 4) for c in b3])
        self.Lfin = np.array([-11 / 6] + [6 * self._fact(i - 1) / self._fact(i + 3) for i in range(1, d2 + 1)]) if r == 3 else None
        self.b2, self.b3 = b2, b3

    @staticmethod
    def _fact(n):
        from math import factorial
        return factorial(n)

    def rows_R1(self, idx, s2, s3, delta, C):
        e = delta + C * self.u[idx]
        R = np.zeros((len(idx), self.nv))
        R[:, :self.n2] = (1 - e * s2)[:, None] * self.B2[idx]
        R[:, self.n2:self.n2 + self.n3] = (1 - e * s3)[:, None] * self.B3[idx]
        sF = np.sign(self.FWc[idx])
        R[:, -1] = (1 - e * sF) * self.FWc[idx]
        return R

    def eval_parts(self, v):
        x = v[:self.n2]; y = v[self.n2:self.n2 + self.n3]; a = v[-1]
        return self.B2 @ x, self.B3 @ y, a * self.FWc

    def matter_rows(self, nm=201):
        """coefficient of a0 (m^0..) and m^4 x coefficient of a2 as linear rows (must be <= 0)."""
        rows = []
        ms = np.linspace(0, 1, nm)
        b2, b3 = self.b2, self.b3
        for m in ms:
            r = np.zeros(self.nv)
            r[:self.n2] = [mom(c, 2) - 2 * m * m * mom(c, 0) for c in b2]
            r[self.n2:self.n2 + self.n3] = [-mom(c, 0) for c in b3]
            rows.append(r)
            r = np.zeros(self.nv)
            m2 = m * m
            r[:self.n2] = [mom(c, 2) * m2 * m2 - 6 * m2 * mom(c, 4) + 6 * mom(c, 6) - 2 * m2 * (m2 * m2 * mom(c, 0) - 6 * m2 * mom(c, 2) + 6 * mom(c, 4)) for c in b2]
            r[self.n2:self.n2 + self.n3] = [-(m2 * m2 * mom(c, 0) - 6 * m2 * mom(c, 2) + 6 * mom(c, 4)) for c in b3]
            r[-1] = -12
            rows.append(r)
        return np.array(rows)

def _rn(R):
    return R / np.maximum(np.abs(R).max(axis=1), 1e-300)[:, None]

def solve(G, delta, C=0.0, u0=1.0, dsc=None, LREF=np.inf, matter=True, maxit=40, verbose=False, vbound=1e5, warm=None):
    """Returns (status, v, info)."""
    if dsc is None:
        dsc = delta
    n2, n3, nv = G.n2, G.n3, G.nv
    # equality: k = 1
    Aeq = np.zeros((1, nv)); Aeq[0, :n2] = -2 * G.mu2_6; Aeq[0, n2:n2 + n3] = 2 * G.mu3_4; Aeq[0, -1] = 4
    beq = [1.0]
    # objective
    c = np.zeros(nv)
    if np.isinf(LREF):
        c[:n2] = 16 * G.Lfin * 1e-3; c[0] += 16
    else:
        c[:n2] = 16 * G.Lfin; c[0] += 16 * (np.log(LREF) - np.log(1.5))
    # fixed rows
    fixed = []
    fixed.append(-np.concatenate([G.Bsc, np.zeros((len(G.bsc), n3 + 1))], axis=1))
    bfix = [-dsc * rho(G.bsc)]
    for h in SIGMA:
        A2, A3, Fc = G.asym[h]
        ok = G.au <= u0
        for s2 in (1, -1):
            for s3 in (1, -1):
                e = delta + C * G.au[ok]
                R = np.zeros((ok.sum(), nv))
                R[:, :n2] = (1 - e * s2)[:, None] * A2[ok]
                R[:, n2:n2 + n3] = (1 - e * s3)[:, None] * A3[ok]
                R[:, -1] = (1 - e * np.sign(Fc[ok])) * Fc[ok]
                fixed.append(-R); bfix.append(np.zeros(ok.sum()))
    if matter:
        M = G.matter_rows()
        fixed.append(M); bfix.append(np.zeros(len(M)))
    Afix = np.concatenate(fixed); bfix = np.concatenate(bfix)
    sc = np.maximum(np.abs(Afix).max(axis=1), 1e-300)
    Afix = Afix / sc[:, None]; bfix = bfix / sc
    # initial cut set: subsample with all 4 patterns
    phys = np.where(G.u <= u0 + 1e-15)[0]
    rng = np.random.default_rng(0)
    init = np.unique(np.concatenate([phys[::max(1, len(phys) // 3000)], phys[G.J[phys] <= 12]]))
    cuts = []
    for s2 in (1, -1):
        for s3 in (1, -1):
            cuts.append(-_rn(G.rows_R1(init, np.full(len(init), s2), np.full(len(init), s3), delta, C)))
    if warm is not None:
        cuts += warm
    bounds = [(-vbound, vbound)] * (nv - 1) + [(0, vbound)]
    for it in range(maxit):
        Aub = np.concatenate([Afix] + cuts); bub = np.concatenate([bfix, np.zeros(sum(len(x) for x in cuts))])
        res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds, method='highs')
        if res.status != 0:
            return res.status, None, {'msg': res.message, 'it': it}
        v = res.x
        B2, B3, FW = G.eval_parts(v)
        FF = B2 + B3 + FW
        N = np.abs(B2) + np.abs(B3) + np.abs(FW)
        g = FF - (delta + C * G.u) * N
        g[G.u > u0 + 1e-15] = np.inf
        viol = np.where(g < -1e-6 * np.maximum(N, 1e-300))[0]
        if verbose:
            print('   it %d obj %.5f  rows %d  violated %d  min g/N %.3e' % (it, res.fun, len(Aub), len(viol),
                  np.min(g[np.isfinite(g)] / np.maximum(N[np.isfinite(g)], 1e-300))), flush=True)
        if len(viol) == 0:
            break
        # add the most violated ones (by relative violation), at most 4000 per iteration, plus local neighbours
        rel = g[viol] / np.maximum(N[viol], 1e-300)
        order = viol[np.argsort(rel)][:4000]
        cuts.append(-G.rows_R1(order, np.sign(B2[order]), np.sign(B3[order]), delta, C) / N[order][:, None])
    else:
        return 9, v, {'msg': 'cutting planes did not converge', 'it': it}
    F = Functional.from_vars(v, G.d2, G.d3, G.r)
    A, B, BA = F.AB()
    return 0, v, {'A': A, 'B': B, 'B_AdS': BA, 'k': F.k(), 'it': it, 'a': v[-1], 'cuts': cuts}

def ratios(G, v):
    """min over physical grid rows of FF/N and the per-row 'envelope threshold' FF/(u N)."""
    B2, B3, FW = G.eval_parts(v)
    FF = B2 + B3 + FW; N = np.abs(B2) + np.abs(B3) + np.abs(FW)
    return FF, N

if __name__ == '__main__':
    d2 = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    d3 = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    G = Grid(d2, d3)
    for delta in [0.0, 0.05, 0.1]:
        st, v, info = solve(G, delta, verbose=True)
        print(delta, st, {k: info[k] for k in info if k != 'cuts'})
