"""NOFORWARD_FUNCTIONAL.md section 2: LP for flat functionals F = int psi2 B2 + psi3 B3 (smeared B_2^{(1)}, B_3^{(1)} only; no forward
terms; these sum rules receive no contact contributions, nf_contacts.py).  Floating point (HiGHS): [NE].

Weights: psi2 = (1-p)^r p sum_{i<=d} x_i p^i,  psi3 = (1-p)^r p sum_{j<=d} y_j p^j  (psi3(0) = 0 is necessary, section 2.1 of the .md).
Heavy action (normal form, u = 1/m^2):  FF = B2 + B3,  B2 = int psi2 (2 - u p^2) D_h(1-2up^2),  B3 = sigma_h u int psi3 D_h;  N = |B2| + |B3|.
Low side (a = 0):  k = -2 mu2_6 + 2 mu3_4 (normalised to 1), A = 16 x_0, B = 16 Lfin, matter coefficients C0, m^4 C2 <= 0.
Constraints: (R1) FF >= delta N on the grid (cutting planes over the sign patterns of B2, B3); (R2) 2 int psi2 J0(bhat p) >= delta rho(bhat);
(R3) bhat -> infinity: -psi2''(0)(1 - delta) >= (1 + delta) u |psi3''(0)|/2 for u in [0, 1] (the 1/bhat^3 coefficient, both channels).
Usage: python3 nf_lp.py r d  -> cache/nf/pareto_r{r}_d{d}.json
"""
import sys, os, json, time
import numpy as np
from scipy.optimize import linprog
from rcommon import pmul, peval, mom, flat_parts_basis, jmax_for, scaling_B2, SIGMA, Functional
from lp import rho

HERE = os.path.dirname(os.path.abspath(__file__))
CDIR = os.path.join(HERE, 'cache', 'nf')

def basis(d, r):
    """(1-p)^r p L_i(2p-1), i = 0..d (shifted Legendre: same span as (1-p)^r p^{i+1}, better conditioned); monomial coefficients."""
    from numpy.polynomial import legendre as Lg, polynomial as Pn
    omp = np.array([1.])
    for _ in range(r):
        omp = pmul(omp, [1., -1.])
    out = []
    for i in range(d + 1):
        li = Lg.leg2poly(np.eye(d + 1)[i])                                # L_i(x) in monomials of x
        sh = np.zeros(1)
        for k, ck in enumerate(li):                                       # L_i(2p - 1)
            sh = Pn.polyadd(sh, ck * Pn.polypow([-1., 2.], k))
        out.append(pmul(pmul(omp, [0., 1.]), sh))
    return out

def basis_values(p, d, r):
    """values of (1-p)^r p L_i(2p-1), i = 0..d, in FACTORED form.  (Evaluating the expanded monomial form near p = 1 loses all
    relative accuracy, which is fatal for +- at u ~ 1 and large J, where P^{(0,8)}_{J-4} ~ J^8 at the backward endpoint.)"""
    from numpy.polynomial import legendre as Lg
    return np.array([(1 - p) ** r * p * Lg.legval(2 * p - 1, np.eye(d + 1)[i]) for i in range(d + 1)])

def parts_factored(h, u, Jmax, d, r):
    from rcommon import nodes, nquad, flat_integrals
    n = nquad(Jmax)
    p, w = nodes(n)
    V = basis_values(p, d, r)
    W = np.vstack([w * V * (2 - u * p * p), w * SIGMA[h] * u * V]).T
    Js, R = flat_integrals(h, u, Jmax, W, p)
    return Js, R[:, :d + 1], R[:, d + 1:]

def u_grid():
    return np.unique(np.concatenate([np.logspace(-4, np.log10(0.25), 50), np.linspace(0.25, 1, 46), np.linspace(0.95, 1, 11)]))

def jmax_nf(u):
    return 1200 if u >= 0.25 else jmax_for(u, 200.0, 2500)

class NFGrid:
    def __init__(self, d, r):
        self.d, self.r = d, r
        self.b = basis(d, r); self.n = d + 1; self.nv = 2 * self.n
        fn = os.path.join(CDIR, 'grid_d%d_r%d.npz' % (d, r))
        if os.path.exists(fn):
            z = np.load(fn); self.B2, self.B3, self.u, self.J, self.h = [z[k] for k in ['B2', 'B3', 'u', 'J', 'h']]
        else:
            L = {k: [] for k in ['B2', 'B3', 'u', 'J', 'h']}
            t0 = time.time()
            for u in u_grid():
                for h in SIGMA:
                    Js, B2, B3 = parts_factored(h, u, jmax_nf(u), d, r)
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
        self.Lfin = np.array([sum(c[i] / (i - 1) for i in range(2, len(c))) + (0 if True else 0) for c in self._lfin_basis()])
        self.pp0 = np.array([2 * c[2] for c in self.b])              # phi''(0)

    def _lfin_basis(self):
        out = []
        for c in self.b:
            c = c.copy(); c[1] -= c[1]                                # int (phi - phi'(0) p)/p^2 : drop the p term
            out.append(c)
        return out

    def rows(self, idx, s2, s3, delta):
        R = np.zeros((len(idx), self.nv))
        R[:, :self.n] = (1 - delta * s2)[:, None] * self.B2[idx]
        R[:, self.n:] = (1 - delta * s3)[:, None] * self.B3[idx]
        return R

    def parts(self, v):
        return self.B2 @ v[:self.n], self.B3 @ v[self.n:]

    def matter_rows(self, nm=201):
        rows = []
        for m in np.linspace(0, 1, nm):
            r = np.zeros(self.nv); m2 = m * m
            r[:self.n] = [mom(c, 2) - 2 * m2 * mom(c, 0) for c in self.b]
            r[self.n:] = [-mom(c, 0) for c in self.b]
            rows.append(r)
            r = np.zeros(self.nv)
            r[:self.n] = [mom(c, 2) * m2 * m2 - 6 * m2 * mom(c, 4) + 6 * mom(c, 6) - 2 * m2 * (m2 * m2 * mom(c, 0) - 6 * m2 * mom(c, 2) + 6 * mom(c, 4)) for c in self.b]
            r[self.n:] = [-(m2 * m2 * mom(c, 0) - 6 * m2 * mom(c, 2) + 6 * mom(c, 4)) for c in self.b]
            rows.append(r)
        return np.array(rows)

def m_proxy_nf(G, delta, BA, um=0.25, saf=4.0):
    """linear sufficient proxy for the Lemma-M certificate at bhat >= BA for ALL u <= um (no forward term):
    (1-d)(-psi2''(0)) - (1+d) um |psi3''(0)|/2 >= saf/BA * [ (4/(3pi))(1+d)(2|psi2'''(1)| + 2||psi2''''||_1 + um(|psi3'''(1)| + ||psi3''''||_1
        + |(p^2psi2)'''(1)| + ||(p^2psi2)''''||_1)) + (2/pi)(1+d) um (6||p (2psi2)'''||_1 + 4||(2psi2)''||_1) ]  (auxiliary |.| variables)."""
    from numpy.polynomial import polynomial as Pn
    n, nv = G.n, G.nv
    x, w = np.polynomial.legendre.leggauss(64); p = (x + 1) / 2; w = w / 2
    der = lambda c, k: Pn.polyder(np.asarray(c, float), k) if len(c) > k else np.zeros(1)
    quants = []
    fam = [(G.b, 0, 2 * (4 / (3 * np.pi)), 3, 1.0), (G.b, n, um * (4 / (3 * np.pi)), 3, 1.0),
           ([np.concatenate([[0, 0], c]) for c in G.b], 0, um * (4 / (3 * np.pi)), 3, 1.0)]
    for basis_, off, wt, _, _ in fam:
        r = np.zeros(nv)
        for i, cc in enumerate(basis_): r[off + i] = Pn.polyval(1.0, der(cc, 3))
        quants.append((r, wt))
        for k in range(len(p)):
            r = np.zeros(nv)
            for i, cc in enumerate(basis_): r[off + i] = Pn.polyval(p[k], der(cc, 4))
            quants.append((r, wt * w[k]))
    for k in range(len(p)):                                   # (2/pi) um (6 ||p (2psi2)'''|| + 4 ||(2psi2)''||)
        r = np.zeros(nv)
        for i, cc in enumerate(G.b): r[i] = 2 * p[k] * Pn.polyval(p[k], der(cc, 3))
        quants.append((r, (2 / np.pi) * um * 6 * w[k]))
        r = np.zeros(nv)
        for i, cc in enumerate(G.b): r[i] = 2 * Pn.polyval(p[k], der(cc, 2))
        quants.append((r, (2 / np.pi) * um * 4 * w[k]))
    naux = len(quants); rows = []
    for j, (r, wt) in enumerate(quants):
        for sg in (1, -1):
            R = np.zeros(nv + naux); R[:nv] = sg * r; R[nv + j] = -1; rows.append(R)
    for sg in (1, -1):
        R = np.zeros(nv + naux)
        R[:n] = (1 - delta) * G.pp0; R[n:nv] = sg * (1 + delta) * um * G.pp0 / 2        # (1-d) psi2'' +- (1+d) um psi3''/2 (<= -proxy)
        for j, (r, wt) in enumerate(quants): R[nv + j] = saf * (1 + delta) * wt / BA
        rows.append(R)
    return np.array(rows), naux

def solve(G, delta, maxit=60, vbound=1e6, verbose=False, mproxy=None):
    n, nv = G.n, G.nv
    Aeq = np.zeros((1, nv)); Aeq[0, :n] = -2 * G.mu2_6; Aeq[0, n:] = 2 * G.mu3_4
    c = np.zeros(nv); c[:n] = 16 * G.Lfin * 1e-3 + 16 * np.array([cc[1] for cc in G.b])
    fixed, bfix = [], []
    fixed.append(-np.hstack([G.Bsc, np.zeros((len(G.bsc), n))])); bfix.append(-delta * rho(G.bsc))
    for uu in np.linspace(0, 1, 41):                                    # (R3)
        for sg in (1, -1):
            r = np.zeros(nv); r[:n] = (1 - delta) * G.pp0; r[n:] = sg * (1 + delta) * uu * G.pp0 / 2
            fixed.append(r[None, :]); bfix.append(np.zeros(1))          # (1-d) psi2''(0) + (1+d) u |psi3''(0)|/2 <= 0
    M = G.matter_rows(); fixed.append(M); bfix.append(np.zeros(len(M)))
    Afix = np.concatenate(fixed); bfix = np.concatenate(bfix)
    sc = np.maximum(np.abs(Afix).max(axis=1), 1e-300); Afix = Afix / sc[:, None]; bfix = bfix / sc
    init = np.unique(np.concatenate([np.arange(0, len(G.u), max(1, len(G.u) // 3000)), np.where(G.J <= 12)[0]]))
    cuts = []
    for s2 in (1, -1):
        for s3 in (1, -1):
            R = -G.rows(init, np.full(len(init), s2), np.full(len(init), s3), delta)
            cuts.append(R / np.maximum(np.abs(R).max(axis=1), 1e-300)[:, None])
    bounds = [(-vbound, vbound)] * nv
    naux = 0
    if mproxy is not None:
        Mr, naux = m_proxy_nf(G, delta, **mproxy)
        pad = lambda R: np.hstack([R, np.zeros((R.shape[0], naux))])
        Afix = np.concatenate([pad(Afix), Mr / np.maximum(np.abs(Mr).max(axis=1), 1e-300)[:, None]]); bfix = np.concatenate([bfix, np.zeros(len(Mr))])
        cuts = [pad(x) for x in cuts]; Aeq = pad(Aeq); c = np.concatenate([c, np.zeros(naux)]); bounds = bounds + [(0, None)] * naux
    for it in range(maxit):
        Aub = np.concatenate([Afix] + cuts); bub = np.concatenate([bfix, np.zeros(sum(len(x) for x in cuts))])
        res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=[1.0], bounds=bounds, method='highs')
        if res.status != 0:
            return res.status, None, {'msg': res.message}
        v = res.x[:nv]
        B2, B3 = G.parts(v); FF = B2 + B3; N = np.abs(B2) + np.abs(B3)
        g = FF - delta * N
        viol = np.where(g < -1e-6 * np.maximum(N, 1e-300))[0]
        if verbose:
            print('   it %d obj %.4f viol %d min g/N %.2e' % (it, res.fun, len(viol), np.min(g / np.maximum(N, 1e-300))), flush=True)
        if len(viol) == 0:
            break
        order = viol[np.argsort(g[viol] / N[viol])][:4000]
        R = -G.rows(order, np.sign(B2[order]), np.sign(B3[order]), delta) / N[order][:, None]
        cuts.append(np.hstack([R, np.zeros((len(R), naux))]) if naux else R)
    else:
        return 9, v, {'msg': 'no convergence'}
    F = to_functional(G, v)
    A, B, BA = F.AB()
    return 0, v, {'A': A, 'B': B, 'B_AdS': BA}

def to_functional(G, v):
    L = max(len(c) for c in G.b)
    p2 = np.zeros(L); p3 = np.zeros(L)
    for i, cc in enumerate(G.b):
        p2[:len(cc)] += v[i] * cc; p3[:len(cc)] += v[G.n + i] * cc
    return Functional(p2, p3, 0.0)

if __name__ == '__main__':
    r, d = int(sys.argv[1]), int(sys.argv[2])
    MP = None
    if len(sys.argv) > 3 and sys.argv[3].startswith('proxy='):
        MP = dict(BA=float(sys.argv[3].split('=')[1])); sys.argv.pop(3)
    deltas = [float(x) for x in sys.argv[3:]] or [0.0, 0.05, 0.1, 0.15, 0.2, 0.3]
    os.makedirs(CDIR, exist_ok=True)
    G = NFGrid(d, r)
    out = []
    for delta in deltas:
        st, v, info = solve(G, delta, mproxy=MP)
        if st != 0:
            print('r=%d d=%d delta=%.2f: status %d %s' % (r, d, delta, st, info.get('msg')), flush=True); continue
        print(('proxy BA=%g ' % MP['BA'] if MP else '') + 'r=%d d=%d delta=%.2f: A=%.4f B=%.4f B_AdS=%.4f' % (r, d, delta, info['A'], info['B'], info['B_AdS']), flush=True)
        out.append(dict(r=r, d=d, delta=delta, v=list(map(float, v)), **info))
        json.dump(out, open(os.path.join(CDIR, 'pareto_r%d_d%d%s.json' % (r, d, '' if MP is None else '_proxy%g' % MP['BA'])), 'w'))
