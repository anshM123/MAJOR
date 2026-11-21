"""Result Q, Stage 7 (RESULT_Q_STAGE7.md §4.3): Stage 5's large-beta positivity re-run with the REALISED forward images of actual weights.

Stage 5 used the exact continuum images (P0 = P_{-1/2}(cosh beta), -D(beta), sigma_h P0).  A weight chi realises
    image(beta) = (sqrt2/pi) int_0^beta W(t) (cosh beta - cosh t)^{-1/2} dt,  W = Phi (1 + e(t)),
so the realised P0, D, C5 images are the Mehler sums of Stage 5 with the extra factor (1 + e(t_i)).  e(t) is computed here in mpmath
(direct quadrature over the weight's own nu grid, 450-bit) on a t grid [0, 30.5] and interpolated (log-linear in |e|, sign kept).
Weights: 'stored' = functionals_q4_deshift.json (EE3 chi_value, chi_deriv) + functionals_q5_deshift.json; 'xp' = Stage-7 extended-precision
weights (cache/stage7_xprec_c4v_d004_N3072.pkl, ..._c4d_..., and the best C5 run).  The scan is stage5_positivity's (same states, beta grid,
normalisation Nrm), forward pieces in 'pure' mode, with  B4 -> 2u^2 P0v, dB4 -> -2u^2 Lam^2 Dd + c4 u^3 P0d, B5 -> sigma u^3 P05.
Output: cache/stage7_s5check.txt
"""
import os, sys, io, time, json, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
import sympy as sp
import stage4_hj as HJ4
import stage5_hj as HJ5
import stage12_forward_regge as FR
import stage5_positivity as S5P

mp.mp.prec = 450
CACHE = os.path.join(HERE, 'cache')
buf = io.StringIO()
TG = np.concatenate([np.linspace(0.05, 6, 25), np.linspace(6.25, 30.5, 98)])


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


class Wc4:
    def __init__(self, chi, dnu, N):
        F = {j: sp.lambdify(HJ4.nu_s, v, 'mpmath') for j, v in HJ4.load()['EE3'].items()}
        self.nu = [(m + mp.mpf(1)/2)*dnu for m in range(N)]
        R = {j: [F[j](mp.mpc(x)) for x in self.nu] for j in F}
        self.A1 = [(R[1][m] + R[-1][m]).real for m in range(N)]; self.A3 = [(R[3][m] + R[-3][m]).real for m in range(N)]
        self.B1 = [(-1j*(R[1][m] - R[-1][m])).real for m in range(N)]; self.B3 = [(-1j*(R[3][m] - R[-3][m])).real for m in range(N)]
        self.chi, self.dnu, self.N = chi, dnu, N

    def __call__(self, t):
        t = mp.mpf(t); ch, sh, ch3, sh3 = mp.cosh(t), mp.sinh(t), mp.cosh(3*t), mp.sinh(3*t)
        tot = mp.mpf(0)
        for m in range(self.N):
            x = self.chi[m]
            if x == 0:
                continue
            c, s = mp.cos(self.nu[m]*t), mp.sin(self.nu[m]*t)
            tot += x*(ch*self.A1[m]*c - sh*self.B1[m]*s + ch3*self.A3[m]*c - sh3*self.B3[m]*s)
        return tot*self.dnu


class Wc5:
    def __init__(self, c, dnu, N):
        F = {j: sp.lambdify(HJ5.nu_s, v, 'mpmath') for j, v in HJ5.load()['C5_pp'].items()}
        self.nu = [(m + mp.mpf(1)/2)*dnu for m in range(N)]
        R = {j: [F[j](mp.mpc(x)) for x in self.nu] for j in F}
        Ns = [R[-2][m] + R[0][m] + R[2][m] for m in range(N)]
        self.hh = {j: [R[j][m]/Ns[m] for m in range(N)] for j in R}
        self.c, self.dnu, self.N = c, dnu, N

    def __call__(self, t):
        t = mp.mpf(t); chs = {j: mp.cosh(j*t) for j in (0, 2, -2)}; shs = {j: mp.sinh(j*t) for j in (0, 2, -2)}
        tot = mp.mpf(0)
        for m in range(self.N):
            cc, ss = mp.cos(self.nu[m]*t), mp.sin(self.nu[m]*t)
            acc = 0
            for j in (0, 2, -2):
                acc += self.hh[j][m]*(cc*chs[j] + 1j*ss*shs[j])
            tot += (self.c[m]*acc).real
        return tot*self.dnu


def e_curve(Wf, Phi):
    out = []
    for t in TG:
        w = Wf(t); ph = Phi(mp.mpf(t))
        out.append(float(w/ph - 1))
    return np.array(out)


def e_interp(ec):
    le = np.log(np.maximum(np.abs(ec), 1e-300)); sg = np.sign(ec)
    def f(t):
        t = np.asarray(t, float)
        v = np.exp(np.interp(t, TG, le, left=le[0], right=le[-1]))
        s = np.interp(t, TG, sg, left=sg[0], right=sg[-1])
        return np.where(t > TG[-1], np.inf, v*np.where(s >= 0, 1, -1))
    return f


def curves():
    fn = os.path.join(CACHE, 'stage7_s5check_curves.pkl')
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    C = {}
    st4 = json.load(open(os.path.join(HERE, 'functionals_q4_deshift.json')))['EE3']
    st5 = json.load(open(os.path.join(HERE, 'functionals_q5_deshift.json')))
    one = lambda t: mp.mpf(1); quad = lambda t: -t*t/2
    d002 = mp.mpf(1)/50; d004 = mp.mpf(1)/25
    t0 = time.time()
    C[('stored', 'c4v')] = e_curve(Wc4([mp.mpf(x) for x in st4['chi_value']], d002, 2048), one)
    C[('stored', 'c4d')] = e_curve(Wc4([mp.mpf(x) for x in st4['chi_deriv']], d002, 2048), quad)
    C[('stored', 'c5')] = e_curve(Wc5([mp.mpc(a, b) for a, b in zip(st5['c_re'], st5['c_im'])], d002, 2048), one)
    pr(f"  stored curves [{time.time() - t0:.0f}s]")
    for piece, Phi in (('c4v', one), ('c4d', quad)):
        d = pickle.load(open(os.path.join(CACHE, f'stage7_xprec_{piece}_d004_N3072.pkl'), 'rb'))
        C[('xp', piece)] = e_curve(Wc4([mp.mpf(x) for x in d['weights']], d004, 3072), Phi)
    c5 = sorted([f for f in os.listdir(CACHE) if f.startswith('stage7_xprec_c5_d004_N2048')])
    best = None
    for f in c5:
        d = pickle.load(open(os.path.join(CACHE, f), 'rb'))
        e30 = [e for t, _, e in d['offgrid'] if t <= 20.2]
        if best is None or max(e30) < best[0]:
            best = (max(e30), f, d)
    pr(f"  C5 xp weight used: {best[1]}")
    C[('xp', 'c5')] = e_curve(Wc5([mp.mpmathify(x) for x in best[2]['weights']], d004, 2048), one)
    C['c5file'] = best[1]
    pr(f"  xp curves [{time.time() - t0:.0f}s]")
    pickle.dump(C, open(fn, 'wb'))
    return C


def scan(Lam, name, betas, us, BQ, E):
    F, c2, c3 = S5P.load_func(name); a, c = F['a'], F['c']
    C4 = {+1: 1.0, -1: 41.0}
    best = {'b>=20': (np.inf, None), 'beta>=2': (np.inf, None), 'beta>=6': (np.inf, None), 'beta>=12': (np.inf, None)}
    first_neg = None
    for beta in betas:
        I20, I22, I3, (P0v, Dd, P0d, P05) = BQ[(Lam, round(beta, 12))]
        bh = Lam*beta
        for h in (+1, -1):
            jmin = 0 if h > 0 else 4
            twoJ1 = Lam*np.sqrt(2*(np.cosh(beta) - 1)/us)
            phys = twoJ1 >= 2*jmin + 1
            B2 = 2*I20 - us*I22; B3 = h*us*I3
            FW = a*(2*us**2*Lam**2*Dd - C4[h]*us**3*P0d + h*us**3*P05) - 2*c*us**2*P0v
            FF = B2 + B3 + FW
            Nrm = np.abs(B2) + np.abs(B3) + a*(2*us**2*Lam**2*abs(Dd) + 42*us**3*abs(P0d) + us**3*abs(P05)) + 2*c*us**2*abs(P0v)
            rat = np.where(phys, FF/Nrm, np.inf)
            i = int(np.argmin(rat))
            if not np.isfinite(rat[i]):
                continue
            J = 0.5*(twoJ1[i] - 1)
            loc = (rat[i], (beta, us[i], h, J, Lam/np.sqrt(us[i]), bh))
            for k, ok in (('b>=20', bh >= 20), ('beta>=2', beta >= 2 and bh >= 20), ('beta>=6', beta >= 6 and bh >= 20),
                          ('beta>=12', beta >= 12 and bh >= 20)):
                if ok and rat[i] < best[k][0]:
                    best[k] = loc
            if first_neg is None and bh >= 20 and rat[i] < 0:
                first_neg = loc
    best['first_neg'] = first_neg if first_neg else (np.inf, None)
    return best


def main():
    t0 = time.time()
    C = curves()
    pr("# stage7_s5check: realised-image check of Stage 5's large-beta positivity")
    for src in ('stored', 'xp'):
        for piece in ('c4v', 'c4d', 'c5'):
            ec = C[(src, piece)]
            pr(f"  {src:6s} {piece}: |W/Phi - 1| at t = 6, 10, 12, 14, 16, 20, 25, 30: " +
               " ".join(f"{abs(np.interp(t, TG, ec)):.1e}" for t in (6, 10, 12, 14, 16, 20, 25, 30)))
    us = np.concatenate([np.geomspace(1e-6, 0.05, 40), np.linspace(0.06, 1.0, 48)])
    for Lam in (100., 300., 1000.):
        betas = np.unique(np.concatenate([np.geomspace(5/Lam, 1.0, 70), np.linspace(1.0, 12.0, 221), np.linspace(12.0, 30.0, 91)]))
        for name in ('Q0-full-robust', 'Q3-full-delta0.2'):
            Fn, c2n, c3n = S5P.load_func(name)
            Tn = [S5P.CosT(c2n), S5P.CosT(np.convolve(c2n, [0, 0, 1])), S5P.CosT(c3n)]
            base = {}
            for beta in betas:
                n = int(min(2.4e5, 800 + 5*Lam*beta)); n = max(800, (n//20)*20)
                t, W = FR.mehler_nodes(beta, n)
                om = Lam*t
                base[beta] = ((W*Tn[0](om)).sum(), (W*Tn[1](om)).sum(), (W*Tn[2](om)).sum(), t, W)
            for src in ('exact', 'stored', 'xp'):
                BQ = {}
                for beta in betas:
                    I20, I22, I3, t, W = base[beta]
                    if src == 'exact':
                        ev = ed = e5 = 0.0
                    else:
                        fv, fd, f5 = (e_interp(C[(src, p)]) for p in ('c4v', 'c4d', 'c5'))
                        ev, ed, e5 = fv(t), fd(t), f5(t)
                    P0v = (W*(1 + ev)).sum(); Dd = 0.5*(W*t**2*(1 + ed)).sum(); P0d = (W*(1 + ed)).sum(); P05 = (W*(1 + e5)).sum()
                    BQ[(Lam, round(beta, 12))] = (I20, I22, I3, (P0v, Dd, P0d, P05))
                best = scan(Lam, name, betas, us, BQ, None)
                pr(f"   Lam = {Lam:5.0f} {name:17s} {src:6s}: " + "; ".join(f"min {k}: {S5P.fmt(v) if v[1] else 'n/a'}" for k, v in best.items()))
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(CACHE, 'stage7_s5check.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
