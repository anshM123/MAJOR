r"""
fn_spin (F_NORM.md, sec. 3.4): the calibration identity in position space, on the spinning-external product family of LIGHTSIDE_DIST
(T' = sigma x T_free, Result Q's swapped functional Psi_k in light-ray polarisation, C x C contour), at the borderline k = 2.

Claim tested (H_EFT(a) at leading order + Lemma A, i.e. the mechanism of the (F) calibration argument):
  for a heavy scalar sigma-exchange W_h (s + t + u channels, dimension Delta_h, external Delta_sigma = 2, d = 3), Regge spin
  n_eff = 0 - 2 < k = 2, so Omega_2[W_h] = S_h + DT[W_h] = 0 (Lemma A; S_h = single-trace (heavy block) residues).  Its EFT (1/m)
  expansion  M_h = -sum_j mu_j (d12^j + d14^j + d13^j),  mu_j = sum_m R_m / delta_m^{j+1},  has its first borderline term at degree 4,
  whose arc depends only on its top Regge coefficient (LIGHTSIDE_DIST (P), [DV]):
       A4 = -(mu_4/2) V[r4],      r4 = (d12 - d14)^4/4
       A5 = -(5/4) mu_5 V[r5],    r5 = (2 - d13)(d12 - d14)^4/4      (degree-5 top Regge term: (5/16)(2 - d13) D^4)
  Prediction:  DT[W_h] = -S_h = A4 (1 + O(1/<delta>)) = A4 + A5 + O(1/<delta>^2),  <delta> ~ mu_4/mu_5 ~ Delta_h^2/O(10).
Numerics: DT[W_h] is evaluated as DT[X4], X4 := W_h - P3 (P3 = degree <= 3 EFT polynomial, n_eff < 2: term-wise value 0 [D; DV in
LIGHTSIDE_DIST sec. 3.2]); X4 = sum_ch sum_m R_m x^4 / (delta_m^4 (x - delta_m)), which removes the large subcritical cancellations.
Exchange residues (s-channel poles delta12 = delta_m = 2 - Delta_h/2 - m):  R_m = ((Delta_h/2 - 1)_m)^2 / (m! (Delta_h - 1/2)_m)
(scalar block, d = 3, Delta_i = 2; overall constant irrelevant: everything is linear in R).
Term-wise residues at delta12 = -n exactly as lsd_run.profile; d13 line Re = 2/3, |Im| <= 8, h = 0.2; production grid of LIGHTSIDE_DIST.

Usage:
  (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 fn_spin.py run ab P1 400)   # and 'ba' in parallel (two processes)  ~6 min each
  python3 fn_spin.py report P1 400                                            -> cache/fn_spin_P1.txt
Read-only imports: lsd_core, lsd_mellin, lsd_pos (and og_free through lsd_core).
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import loggamma, digamma, gammaln

CACHE = os.path.join(HERE, 'cache')
GRID = dict(s0=1e-9, m=16, q=1.8)
DHS = (13.0, 21.0, 29.0, 41.0, 57.0)


def m_nodes(Dh, M0=6000, nq=200):
    """nodes/weights for sums over m >= 0 of smooth summands times R_m: direct m < M0, tail int_{M0-1/2}^oo (Gauss-Legendre, mapped)
    + midpoint Euler-Maclaurin derivative term.  Returns (m, w*R_m normalised, delta_m, log-norm)."""
    md = np.arange(M0, dtype=float)
    x, w = np.polynomial.legendre.leggauss(nq)
    uu = 0.5*(x + 1); wu = 0.5*w
    x0 = M0 - 0.5; L = 4.0*M0
    mi = x0 + L*uu/(1 - uu); wi = wu*L/(1 - uu)**2
    st = x0 + np.array([-2.0, -1.0, 1.0, 2.0]); wst = np.array([1.0, -8.0, 8.0, -1.0])/12.0/24.0
    m = np.concatenate([md, mi, st]); wm = np.concatenate([np.ones(M0), wi, wst])
    a = Dh/2 - 1; b = Dh - 0.5
    lR = 2*(gammaln(a + m) - gammaln(a)) - gammaln(m + 1) - (gammaln(b + m) - gammaln(b))
    L0 = lR.max()
    return m, wm*np.exp(lR - L0), 2 - Dh/2 - m, L0


class Exch:
    def __init__(self, Dh):
        self.Dh = Dh
        self.m, self.wR, self.dm, self.L0 = m_nodes(Dh)
        self.mu = np.array([np.sum(self.wR/self.dm**(j + 1)) for j in range(8)])

    def X(self, x):
        """X4(x) = sum_m R_m x^4/(delta_m^4 (x - delta_m)) and dX/dx, for an array x."""
        x = np.asarray(x, complex)[:, None]
        dm = self.dm[None, :]
        q = 1.0/(x - dm)
        base = x**4/dm**4
        Xv = np.sum(self.wR[None, :]*base*q, axis=1)
        dX = np.sum(self.wR[None, :]*(4*x**3/dm**4*q - base*q*q), axis=1)
        return Xv, dX


def run(order, tag, N=400, Y13=8.0, h13=0.2, c2=2/3):
    import lsd_core as C
    import lsd_pos as P
    a, b, p2, p4 = P.config(tag)
    aa, bb = (a, b) if order == 'ab' else (b, a)
    t0 = time.time()
    S = C.Setup(2, aa, bb, p2, p4, **GRID)
    print(f'setup {len(S.t)} nodes [{time.time() - t0:.0f}s]', flush=True)
    EX = {Dh: Exch(Dh) for Dh in DHS}
    nf = np.arange(N + 1).astype(float)
    dL = -nf.astype(complex)
    pieces = ['r4', 'r5'] + [f'X{Dh:g}' for Dh in DHS]
    out = {(p, c): np.zeros(N + 1, complex) for p in pieces for c in S.H}
    y13 = np.arange(-Y13, Y13 + h13/2, h13)
    for i, yy in enumerate(y13):
        d13 = c2 + 1j*yy
        w13 = np.exp(2*loggamma(d13))*h13/(2*np.pi)
        Kd = S.K(dL, d13, deriv=True)
        gL = np.exp(2*loggamma(2 - d13 + nf) - 2*loggamma(nf + 1)); psi14 = digamma(2 - d13 + nf); psin = digamma(nf + 1)
        d14 = 2 - dL - d13
        D = dL - d14
        Ms = {'r4': (D**4/4, 2*D**3), 'r5': ((2 - d13)*D**4/4, (2 - d13)*2*D**3)}
        for Dh, E in EX.items():
            X12, dX12 = E.X(dL)
            X14, dX14 = E.X(d14)
            X13, _ = E.X(np.array([d13]))
            Ms[f'X{Dh:g}'] = (X12 + X14 + X13[0], dX12 - dX14)
        for p, (Mv, dMv) in Ms.items():
            for c, (K, dK) in Kd.items():
                out[(p, c)] += w13*gL*((-2*psi14*Mv + dMv)*K + Mv*dK + 2*psin*Mv*K)
        if i % 10 == 0:
            print(f'  d13 node {i}/{len(y13)} [{time.time() - t0:.0f}s]', flush=True)
    res = {'order': order, 'tag': tag, 'N': N, 'grid': GRID, 'res': out,
           'mu': {Dh: E.mu for Dh, E in EX.items()}, 'L0': {Dh: E.L0 for Dh, E in EX.items()}}
    with open(os.path.join(CACHE, f'fn_spin_{tag}_{order}_N{N}.pkl'), 'wb') as fh:
        pickle.dump(res, fh)
    print(f'done [{time.time() - t0:.0f}s]', flush=True)


def tail(r, n0):
    n = np.arange(n0, len(r)); y = r[n0:]
    pp = -np.polyfit(np.log(n), np.log(np.abs(y) + 1e-300), 1)[0]
    A = np.sum(y*n**pp)/len(n)
    if pp <= 1.05:
        return pp, np.nan
    return pp, A*(len(r) - 0.5)**(1 - pp)/(pp - 1)


def report(tag, N=400, NU=160, NT=80):
    """NU: residues used (n <= NU; the production grid resolves K(-n) to n ~ 170, cf. LIGHTSIDE_DIST's N = 150); tail fit on [NT, NU]."""
    L = []
    pr = lambda s: (print(s, flush=True), L.append(s))
    R = {o: pickle.load(open(os.path.join(CACHE, f'fn_spin_{tag}_{o}_N{N}.pkl'), 'rb')) for o in ('ab', 'ba')}
    cls = sorted(set(c for (p, c) in R['ab']['res']))
    def symarr(p):
        r = sum(R['ab']['res'][(p, c)] for c in cls) + sum(R['ba']['res'][(p, c)] for c in cls)     # k = 2: sym = ab + ba
        return r[:NU + 1]
    pr(f'fn_spin report, tag {tag}, k = 2; residues n <= {NU} (grid-resolved range), tail fit on n in [{NT}, {NU}]; sym = (ab) + (ba),'
       f' six free classes')
    V = {}
    for p in ('r4', 'r5'):
        r = symarr(p); pp, tl = tail(r, NT)
        V[p] = np.sum(r) + tl
        pr(f'  {p}: sym N=100 {np.sum(r[:101]).real:+.7e}, N=150 {np.sum(r[:151]).real:+.7e}, N={NU} {np.sum(r).real:+.7e}; '
           f'tail {tl.real:+.2e} (p = {pp:.2f}) -> V[{p}] = {V[p].real:+.7e}')
    pr('  heavy scalar sigma-exchange, X4 = W_h - P3 (term-wise value = DT[W_h] = -S_h by Lemma A); units: R_m normalised to max_m R_m = 1')
    pr('  Delta_h  <delta>     DT[X4] (sym+tail)    A4 = -(mu4/2)V[r4]   A4+A5               DT/A4-1     DT/(A4+A5)-1  | tail/DT, p; DT(N=120)/DT(N=160)-1')
    rows = []
    for Dh in DHS:
        mu = R['ab']['mu'][Dh]
        r = symarr(f'X{Dh:g}'); pp, tl = tail(r, NT)
        tl = tl if np.isfinite(tl) else 0.0
        DT = np.sum(r) + tl
        A4 = -(mu[4]/2)*V['r4']; A5 = -(5/4)*mu[5]*V['r5']
        rows.append((Dh, (mu[4]/mu[5]).real, DT, A4, A5))
        pr(f'  {Dh:6.1f} {(mu[4]/mu[5]).real:9.3f}   {DT.real:+.8e}   {A4.real:+.8e}   {(A4 + A5).real:+.8e}   '
           f'{(DT/A4 - 1).real:+.3e}  {(DT/(A4 + A5) - 1).real:+.3e}  | {abs(tl)/abs(DT):.1e}, {pp:.2f}; {(np.sum(r[:121])/np.sum(r) - 1).real:+.1e}')
    d = np.array([x[1] for x in rows]); e4 = np.array([(x[2]/x[3] - 1).real for x in rows]); e45 = np.array([(x[2]/(x[3] + x[4]) - 1).real for x in rows])
    pr(f'  scaling: <delta> x (DT/A4 - 1) = ' + ', '.join(f'{v:+.3f}' for v in d*e4) + ';   <delta>^2 x (DT/(A4+A5) - 1) = ' +
       ', '.join(f'{v:+.2f}' for v in d**2*e45))
    with open(os.path.join(CACHE, f'fn_spin_{tag}.txt'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')


def run_sh(order, tag, M=160, Y13=8.0, h13=0.2, c2=2/3, dhs=(13.0, 21.0)):
    """single-trace (heavy block) residues S_h: at delta12 = delta_m = 2 - Delta_h/2 - m, Res = Gamma(delta_m)^2 Gamma(2 - delta_m - d13)^2
    R_m K(delta_m, d13) (s-channel pole of W_h; X4 has the same poles and residues).  m = 0..M (direct nodes of m_nodes)."""
    import lsd_core as C
    import lsd_pos as P
    a, b, p2, p4 = P.config(tag)
    aa, bb = (a, b) if order == 'ab' else (b, a)
    t0 = time.time()
    S = C.Setup(2, aa, bb, p2, p4, **GRID)
    EX = {Dh: Exch(Dh) for Dh in dhs}
    out = {(Dh, c): np.zeros(M + 1, complex) for Dh in dhs for c in S.H}
    y13 = np.arange(-Y13, Y13 + h13/2, h13)
    for i, yy in enumerate(y13):
        d13 = c2 + 1j*yy
        w13 = np.exp(2*loggamma(d13))*h13/(2*np.pi)
        for Dh, E in EX.items():
            dm = E.dm[:M + 1].astype(complex)
            Kd = S.K(dm, d13)
            g = np.exp(2*loggamma(dm) + 2*loggamma(2 - dm - d13))*E.wR[:M + 1]
            for c, K in Kd.items():
                out[(Dh, c)] += w13*g*K
        if i % 20 == 0:
            print(f'  d13 node {i}/{len(y13)} [{time.time() - t0:.0f}s]', flush=True)
    with open(os.path.join(CACHE, f'fn_spin_sh_{tag}_{order}.pkl'), 'wb') as fh:
        pickle.dump({'res': out, 'dhs': dhs, 'M': M}, fh)
    print(f'done [{time.time() - t0:.0f}s]', flush=True)


def report_sh(tag, N=600, NU=160, NT=80):
    R = {o: pickle.load(open(os.path.join(CACHE, f'fn_spin_{tag}_{o}_N{N}.pkl'), 'rb')) for o in ('ab', 'ba')}
    Q = {o: pickle.load(open(os.path.join(CACHE, f'fn_spin_sh_{tag}_{o}.pkl'), 'rb')) for o in ('ab', 'ba')}
    cls = sorted(set(c for (p, c) in R['ab']['res']))
    L = []
    pr = lambda s: (print(s, flush=True), L.append(s))
    pr(f'fn_spin Lemma-A check, tag {tag}, k = 2: S_h (single-trace sum, m <= {Q["ab"]["M"]}) + DT[X4] (n <= {NU}) per order and sym')
    for Dh in Q['ab']['dhs']:
        for o in ('ab', 'ba'):
            sh = sum(Q[o]['res'][(Dh, c)] for c in cls)
            dt = sum(R[o]['res'][(f'X{Dh:g}', c)][:NU + 1] for c in cls)
            ps, ts = tail(sh, NT); pd, td = tail(dt, NT)
            ts = ts if np.isfinite(ts) else 0; td = td if np.isfinite(td) else 0
            Sv = np.sum(sh) + ts; Dv = np.sum(dt) + td
            pr(f'  Delta_h = {Dh}, ({o}): S_h = {Sv.real:+.8e} (tail {abs(ts)/abs(Sv):.1e}, p {ps:.2f}), DT = {Dv.real:+.8e} '
               f'(tail {abs(td)/abs(Dv):.1e}), S_h + DT = {(Sv + Dv).real:+.3e}  (rel. {abs(Sv + Dv)/max(abs(Sv), abs(Dv)):.1e})')
    with open(os.path.join(CACHE, f'fn_spin_sh_{tag}.txt'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(sys.argv[2], sys.argv[3], N=int(sys.argv[4]) if len(sys.argv) > 4 else 400)
    elif sys.argv[1] == 'report':
        report(sys.argv[2], N=int(sys.argv[3]) if len(sys.argv) > 3 else 400)
    elif sys.argv[1] == 'runsh':
        run_sh(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == 'reportsh':
        report_sh(sys.argv[2])
    elif sys.argv[1] == 'mu':
        for Dh in DHS:
            E = Exch(Dh)
            print(Dh, 'mu4/mu5 =', E.mu[4]/E.mu[5], ' mu3/mu4 =', E.mu[3]/E.mu[4], ' mu5/mu6 =', E.mu[5]/E.mu[6])
