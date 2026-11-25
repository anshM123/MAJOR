"""LF1/LF2: flat LP for a low-nu cubic functional.  Units M' = a*Lam (the support edge): p' in [0,1], heavy states u' = (a Lam/m)^2 in (0, a^2].
By the dilation identity (LOWNU_SCOPE section 2) this is exactly the problem "weights supported on nu <= a Lam, positivity on every m >= Lam".
Same weight class and constraints as robust/nf_lp.py (psi_k = (1-p)^r p * Legendre_d, (R1) FF >= delta N by cutting planes on a (u', J, h) grid,
(R2) scaling limit, (R3) 1/bhat^3 coefficient for u' <= a^2, optional matter rows for light scalars/spin-2 with m_l <= 1/a, i.e. below Lam),
objective = 16 L_AdS(Lam') = 16[x0 log(2 Lam'/3) + Lfin] at normalisation k = 1 (Lam' = a Lam).   Floating point (HiGHS): [NE].
Usage: python3 ln2_lp.py a r d delta [matter=0/1] [Lam=998.5]"""
import sys, os, time, json
import numpy as np
from scipy.optimize import linprog
from ln_common import guard, CACHE
import ln_common  # noqa  (sets path)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'robust'))
from rcommon import mom, jmax_for, scaling_B2, SIGMA  # read-only
from lp import rho                                    # read-only
from nf_lp import basis, parts_factored, to_functional  # read-only (no writes at import)

def ugrid(a):
    a2 = a * a
    return np.unique(np.concatenate([np.logspace(-5, np.log10(a2), 36), np.linspace(0.5 * a2, a2, 11)]))

class Grid:
    def __init__(self, a, d, r):
        self.a, self.d, self.r = a, d, r
        self.b = basis(d, r); self.n = d + 1; self.nv = 2 * self.n
        fn = guard(os.path.join(CACHE, 'grid_a%.4f_d%d_r%d.npz' % (a, d, r)))
        if os.path.exists(fn):
            z = np.load(fn); self.B2, self.B3, self.u, self.J, self.h = [z[k] for k in ['B2', 'B3', 'u', 'J', 'h']]
        else:
            L = {k: [] for k in ['B2', 'B3', 'u', 'J', 'h']}; t0 = time.time()
            for u in ugrid(a):
                for h in SIGMA:
                    Js, B2, B3 = parts_factored(h, u, jmax_for(u, 200.0, 2500), d, r)
                    if h == '++':
                        sel = Js % 2 == 0; Js, B2, B3 = Js[sel], B2[sel], B3[sel]
                    L['B2'].append(B2); L['B3'].append(B3); L['u'].append(np.full(len(Js), u)); L['J'].append(Js)
                    L['h'].append(np.full(len(Js), SIGMA[h]))
            for k in L: setattr(self, k, np.concatenate(L[k]))
            np.savez(fn, B2=self.B2, B3=self.B3, u=self.u, J=self.J, h=self.h)
            print('  grid built: %d rows (%.0fs)' % (len(self.u), time.time() - t0), flush=True)
        self.bsc = np.unique(np.concatenate([np.linspace(0, 60, 1201), np.linspace(60, 400, 1700)]))
        self.Bsc = scaling_B2(self.b, self.bsc)
        self.mu2_6 = np.array([mom(c, 6) for c in self.b]); self.mu3_4 = np.array([mom(c, 4) for c in self.b])
        self.x0 = np.array([c[1] for c in self.b])
        self.Lfin = np.array([sum(c[i] / (i - 1) for i in range(2, len(c))) for c in self.b])
        self.pp0 = np.array([2 * c[2] for c in self.b])
    def rows(self, idx, s2, s3, delta):
        R = np.zeros((len(idx), self.nv))
        R[:, :self.n] = (1 - delta * s2)[:, None] * self.B2[idx]
        R[:, self.n:] = (1 - delta * s3)[:, None] * self.B3[idx]
        return R
    def parts(self, v):
        return self.B2 @ v[:self.n], self.B3 @ v[self.n:]
    def matter_rows(self, mmax, nm=401):
        rows = []
        for m in np.linspace(0, mmax, nm):
            r = np.zeros(self.nv); m2 = m * m
            r[:self.n] = [mom(c, 2) - 2 * m2 * mom(c, 0) for c in self.b]; r[self.n:] = [-mom(c, 0) for c in self.b]; rows.append(r)
            r = np.zeros(self.nv)
            r[:self.n] = [mom(c, 2) * m2 * m2 - 6 * m2 * mom(c, 4) + 6 * mom(c, 6) - 2 * m2 * (m2 * m2 * mom(c, 0) - 6 * m2 * mom(c, 2) + 6 * mom(c, 4)) for c in self.b]
            r[self.n:] = [-(m2 * m2 * mom(c, 0) - 6 * m2 * mom(c, 2) + 6 * mom(c, 4)) for c in self.b]; rows.append(r)
        return np.array(rows)

def solve(G, delta, Lamp, matter, maxit=80, vbound=1e7):
    n, nv = G.n, G.nv; a2 = G.a ** 2
    Aeq = np.zeros((1, nv)); Aeq[0, :n] = -2 * G.mu2_6; Aeq[0, n:] = 2 * G.mu3_4
    from scipy.integrate import quad
    e2 = 9 / (4 * Lamp ** 2)
    LAb = np.array([quad(lambda p: np.polyval(cb[::-1], p) / (p * p + e2), 0, 1, points=[1.5 / Lamp, 0.05], limit=400, epsabs=0, epsrel=1e-12)[0] for cb in G.b])
    c = np.zeros(nv); c[:n] = 16 * LAb
    fixed, bfix = [-np.hstack([G.Bsc, np.zeros((len(G.bsc), n))])], [-delta * rho(G.bsc)]
    for uu in np.linspace(0, a2, 21):
        for sg in (1, -1):
            r = np.zeros(nv); r[:n] = (1 - delta) * G.pp0; r[n:] = sg * (1 + delta) * uu * G.pp0 / 2
            fixed.append(r[None, :]); bfix.append(np.zeros(1))
    if matter:
        M = G.matter_rows(1.0 / G.a); fixed.append(M); bfix.append(np.zeros(len(M)))
    Afix = np.concatenate(fixed); bfix = np.concatenate(bfix)
    sc = np.maximum(np.abs(Afix).max(axis=1), 1e-300); Afix = Afix / sc[:, None]; bfix = bfix / sc
    init = np.unique(np.concatenate([np.arange(0, len(G.u), max(1, len(G.u) // 3000)), np.where(G.J <= 12)[0]]))
    cuts = []
    for s2 in (1, -1):
        for s3 in (1, -1):
            R = -G.rows(init, np.full(len(init), s2), np.full(len(init), s3), delta)
            cuts.append(R / np.maximum(np.abs(R).max(axis=1), 1e-300)[:, None])
    for it in range(maxit):
        Aub = np.concatenate([Afix] + cuts); bub = np.concatenate([bfix, np.zeros(sum(len(x) for x in cuts))])
        res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=[1.0], bounds=[(-vbound, vbound)] * nv, method='highs')
        if res.status != 0:
            return res.status, None, {'msg': res.message}
        v = res.x
        B2, B3 = G.parts(v); FF = B2 + B3; N = np.abs(B2) + np.abs(B3); g = FF - delta * N
        viol = np.where(g < -1e-7 * np.maximum(N, 1e-300))[0]
        if len(viol) == 0:
            break
        order = viol[np.argsort(g[viol] / N[viol])][:4000]
        cuts.append(-G.rows(order, np.sign(B2[order]), np.sign(B3[order]), delta) / N[order][:, None])
    else:
        return 9, v, {'msg': 'no convergence'}
    F = to_functional(G, v)
    x0 = F.psi2[1]; Lfin = F.Lfin()
    return 0, v, {'obj16L': float(res.fun), 'x0': float(x0), 'Lfin': float(Lfin), 'A': 16 * x0, 'B': 16 * Lfin,
                  'minFFN_grid': float(np.min(FF / np.maximum(N, 1e-300))), 'it': it}

if __name__ == '__main__':
    a, r, d, delta = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
    matter = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    Lam = float(sys.argv[6]) if len(sys.argv) > 6 else 998.5
    G = Grid(a, d, r)
    st, v, info = solve(G, delta, a * Lam, matter)
    if st != 0:
        print('a=%.3f r=%d d=%d delta=%.2f matter=%d: status %d %s' % (a, r, d, delta, matter, st, info.get('msg'))); sys.exit()
    # tau from the exact LP weight with 9/(4 Lam'^2) regulator (k = 1): tau = 240 sqrt(16 L(Lam'))/Lam'^4
    from scipy.integrate import quad
    F = to_functional(G, v); Lp = a * Lam; e2 = 9 / (4 * Lp ** 2)
    pv = lambda p: np.polyval(F.psi2[::-1], p)
    LA = quad(lambda p: pv(p) / (p * p + e2), 0, 1, points=[1.5 / Lp, 0.05], limit=400)[0]
    tau = 240 * np.sqrt(16 * LA) / Lp ** 4
    info.update(a=a, r=r, d=d, delta=delta, matter=matter, Lam=Lam, L_A=LA, tau=float(tau), v=list(map(float, v)))
    print('a=%.3f r=%d d=%d delta=%.2f matter=%d: 16L=%.2f (A=%.1f B=%.1f) L_A=%.2f tau=%.4e  minFF/N(grid)=%.4f it=%d' %
          (a, r, d, delta, matter, info['obj16L'], info['A'], info['B'], LA, tau, info['minFFN_grid'], info['it']), flush=True)
    fn = guard(os.path.join(CACHE, 'lp_a%.4f_r%d_d%d_dl%.2f_m%d.json' % (a, r, d, delta, matter)))
    json.dump(info, open(fn, 'w'))
