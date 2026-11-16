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

    def add_u(self, us, bmax=200.0, Jcap=2500):
        """append physical rows at extra u values (e.g. dense near the endpoint u0 of an envelope-targeted LP)."""
        b2 = basis2(self.d2, self.r); b3 = basis3(self.d3, self.r)
        add = {k: [getattr(self, k)] for k in ['B2', 'B3', 'FWc', 'u', 'J', 'h']}
        for u in us:
            Jm = jmax_for(u, bmax, Jcap)
            for h in ['++', '+-']:
                Js, B2, B3 = flat_parts_basis(h, u, Jm, b2, b3)
                if h == '++':
                    sel = Js % 2 == 0; Js, B2, B3 = Js[sel], B2[sel], B3[sel]
                add['B2'].append(B2); add['B3'].append(B3); add['FWc'].append(phi_coef(h, Js, u))
                add['u'].append(np.full(len(Js), u)); add['J'].append(Js); add['h'].append(np.full(len(Js), SIGMA[h]))
        for k, v in add.items():
            setattr(self, k, np.concatenate(v))

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

def LAdS_rows(G, Lams, n=4000):
    """16 int_0^1 phi(p)/(p^2 + 9/(4 Lams^2)) dp for each psi2 basis function (exact AdS light-side log, ADS4_LIGHT.md)."""
    from rcommon import gl01
    # substitution p = a sinh(t)-type clustering is unnecessary: phi ~ p near 0, integrand bounded by 1/(2a); use graded GL panels
    a = 1.5 / Lams
    edges = np.unique(np.concatenate([[0.0], a * np.logspace(-3, np.log10(1 / a), 60), [1.0]]))
    edges = edges[edges <= 1.0]
    x, w = np.polynomial.legendre.leggauss(40)
    tot = np.zeros(G.n2)
    for lo, hi in zip(edges[:-1], edges[1:]):
        p = lo + (x + 1) * (hi - lo) / 2; ww = w * (hi - lo) / 2
        for i, c in enumerate(G.b2):
            tot[i] += np.sum(ww * peval(c, p) / (p * p + a * a))
    return 16 * tot

def solve(G, delta, C=0.0, u0=1.0, dsc=None, LREF=np.inf, matter=True, maxit=40, verbose=False, vbound=1e5, warm=None, LAMS=None, mreg=None):
    """Returns (status, v, info).
    mreg = dict(BA, u1, Y, amin, safety): adds a linear sufficient proxy for the Lemma-M certificate at bhat >= BA (section 4 of the .md):
      (1-delta) gamma0 BA/(1+ds)^3 >= safety (4/(3pi)) [2(1+delta)(|psi2'''(1)| + ||psi2''''||_1) + u1(1+delta)(|psi3'''(1)| + ||psi3''''||_1
                                       + |(p^2psi2)'''(1)| + ||(p^2psi2)''''||_1)]/(1-ds)^4 + 1.2 (1+delta)^2 Y^2/(2(1-delta) amin (1-ds)^2),
      gamma0 = -psi2''(0), |psi3(0)| <= Y, a >= amin  (L1 norms by 64-point Gauss-Legendre, absolute values via auxiliary variables)."""
    if dsc is None:
        dsc = delta
    n2, n3, nv = G.n2, G.n3, G.nv
    # equality: k = 1
    Aeq = np.zeros((1, nv)); Aeq[0, :n2] = -2 * G.mu2_6; Aeq[0, n2:n2 + n3] = 2 * G.mu3_4; Aeq[0, -1] = 4
    beq = [1.0]
    # objective
    c = np.zeros(nv)
    if LAMS is not None:
        c[:n2] = LAdS_rows(G, LAMS)
    elif np.isinf(LREF):
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
    naux = 0
    if mreg is not None:
        Mrows, Mb, naux, bounds = m_proxy(G, delta, mreg, bounds)
    padc = lambda R: np.hstack([R, np.zeros((R.shape[0], naux))]) if naux else R
    Afix = padc(Afix); cuts = [padc(x) for x in cuts]; Aeq = padc(Aeq); c = np.concatenate([c, np.zeros(naux)])
    if mreg is not None:
        Afix = np.concatenate([Afix, Mrows]); bfix = np.concatenate([bfix, Mb])
    for it in range(maxit):
        Aub = np.concatenate([Afix] + cuts); bub = np.concatenate([bfix, np.zeros(sum(len(x) for x in cuts))])
        res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds, method='highs')
        if res.status != 0:
            return res.status, None, {'msg': res.message, 'it': it}
        v = res.x[:nv]
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
        cuts.append(padc(-G.rows_R1(order, np.sign(B2[order]), np.sign(B3[order]), delta, C) / N[order][:, None]))
    else:
        return 9, v, {'msg': 'cutting planes did not converge', 'it': it}
    F = Functional.from_vars(v, G.d2, G.d3, G.r)
    A, B, BA = F.AB()
    return 0, v, {'A': A, 'B': B, 'B_AdS': BA, 'k': F.k(), 'it': it, 'a': v[-1], 'cuts': cuts}

def m_proxy(G, delta, mreg, bounds):
    from numpy.polynomial import polynomial as Pn
    BA, u1, Y, amin, saf = mreg['BA'], mreg['u1'], mreg['Y'], mreg['amin'], mreg.get('safety', 1.4)
    n2, n3, nv = G.n2, G.n3, G.nv
    x, w = np.polynomial.legendre.leggauss(64); p = (x + 1) / 2; w = w / 2
    quants = []   # (row over v, weight in the proxy)
    def der(c, k): return Pn.polyder(np.asarray(c, float), k) if len(c) > k else np.zeros(1)
    fam = [('psi2', G.b2, 0, 2 * (1 + delta)), ('psi3', G.b3, n2, u1 * (1 + delta)),
           ('p2psi2', [np.concatenate([[0, 0], c]) for c in G.b2], 0, u1 * (1 + delta))]
    for nm, basis, off, wt in fam:
        r = np.zeros(nv)
        for i, cc in enumerate(basis): r[off + i] = Pn.polyval(1.0, der(cc, 3))
        quants.append((r, wt))
        for k in range(len(p)):
            r = np.zeros(nv)
            for i, cc in enumerate(basis): r[off + i] = Pn.polyval(p[k], der(cc, 4))
            quants.append((r, wt * w[k]))
    naux = len(quants)
    rows, bb = [], []
    for j, (r, wt) in enumerate(quants):
        for sg in (1, -1):
            R = np.zeros(nv + naux); R[:nv] = sg * r; R[nv + j] = -1; rows.append(R); bb.append(0.0)
    ds = 8 * np.sqrt(u1) / BA
    R = np.zeros(nv + naux)
    g0 = np.zeros(nv); g0[:n2] = [-2 * c[2] for c in G.b2]          # -psi2''(0)
    R[:nv] = -(1 - delta) * g0 * BA / (1 + ds) ** 3
    for j, (r, wt) in enumerate(quants): R[nv + j] = saf * (4 / (3 * np.pi)) * wt / (1 - ds) ** 4
    Ka = 1.2 * (1 + delta) ** 2 * Y ** 2 / (2 * (1 - delta) * amin * (1 - ds) ** 2)
    rows.append(R); bb.append(-Ka)
    bounds = list(bounds)
    bounds[n2] = (-Y, Y)                       # y_0 = psi3(0)
    bounds[nv - 1] = (amin, bounds[nv - 1][1])
    bounds += [(0, None)] * naux
    return np.array(rows), np.array(bb), naux, bounds

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
