"""Result Q, Stage 5 (RESULT_Q_STAGE5.md §4): heavy positivity of the delta = 0.05 (Q0-full-robust) and delta = 0.2 (Q3-full-delta0.2)
functionals in Stage 12's leading-Regge model, extended to LARGE J (beta <= 30, i.e. J up to ~1e8 at m = Lam), all u in [1e-6, 1],
with BOTH forward pieces de-shifted (C4: stage 4; C5: stage 5, spinning) -> pure images
    B4|_0 -> 2u^2 P0,  dB4|_0 -> -2u^2 Lam^2 D + c4_h u^3 P0,  B5|_0 -> sigma_h u^3 P0,
and, for comparison, with the RAW spinning C5 (exact image, stage5_hj):
    B5|_0 -> sigma_h u^3 H5,  H5 = k1 d_mu P_mu|_{3/2} + k2 P_{3/2} + k3 P_{-1/2},
    k1 = -192491/1028120, k2 = 1252043/1028120, k3 = -223923/1028120   (H5(0) = 1; ++ image; +- = -(++)).
Smeared pieces exactly as stage12_forward_regge (b): int psi2 (2 - u p^2) P_{-1/2+i Lam p} + sigma u int psi3 P_{-1/2+i Lam p}, by Mehler-
Dirichlet with EXACT polynomial cos transforms (closed-form boundary sums for omega >= 80, Gauss-Legendre below) so that Lam beta up to 3e4
is resolved.  Physical states: 2J + 1 = Lam sqrt(2(cosh beta - 1)/u) >= 2 J_min + 1 (J_min = 0 (++), 4 (+-)).  Model validity (Stage 12
§1.5): bhat = Lam beta >= 20; bhat in [5, 20) reported separately (outside the model's validity; near-gap region is Stage 3's).
Normalisation as stage12/stage4_c5: Nrm = |B2| + |B3| + a(2u^2 Lam^2 D + 42 u^3 P0 + u^3 |H5|) + 2c u^2 P0.
Output: cache/stage5_positivity.txt
"""
import os, sys, io, json, time
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
import stage12_forward_regge as FR

CACHE = os.path.join(HERE, 'cache')
K1, K2, K3 = -192491/1028120, 1252043/1028120, -223923/1028120
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def load_func(name):
    for fn in ('functionals_q.json', 'functionals_q3.json'):
        D = json.load(open(os.path.join(HERE, fn)))
        if name in D:
            F = D[name]
            return F, np.array(F['psi2'], float), np.array(F['psi3'], float) if len(F['psi3']) else np.zeros(1)
    raise KeyError(name)


class CosT:
    """exact int_0^1 q(p) cos(w p) dp for a polynomial q (coefficients low -> high)."""
    def __init__(self, coef, wc=80.0):
        self.q = np.polynomial.Polynomial(coef)
        self.d0 = []; self.d1 = []
        qq = self.q
        for k in range(len(coef) + 1):
            self.d0.append(qq(0.0)); self.d1.append(qq(1.0)); qq = qq.deriv()
        self.wc = wc
        p, wp = FR.composite_gl(0., 1., 40, 20)
        self.p, self.wq = p, wp*self.q(p)

    def __call__(self, w):
        w = np.asarray(w, float); out = np.empty_like(w)
        lo = w < self.wc
        if np.any(lo):
            out[lo] = np.cos(np.outer(w[lo], self.p)) @ self.wq
        hi = ~lo
        if np.any(hi):
            x = w[hi]; e = np.exp(1j*x); acc = 0
            iw = 1j*x
            for k in range(len(self.d0)):
                acc = acc + (-1)**k*(self.d1[k]*e - self.d0[k])/iw**(k + 1)
            out[hi] = acc.real
        return out


def beta_quantities(beta, Lam, T, n=None):
    if n is None:
        n = int(min(2.4e5, 800 + 5*Lam*beta))
        n = max(800, (n//20)*20)
    t, W = FR.mehler_nodes(beta, n)
    om = Lam*t
    I20 = (W*T[0](om)).sum(); I22 = (W*T[1](om)).sum(); I3 = (W*T[2](om)).sum()
    P0 = W.sum(); D = 0.5*(W*t**2).sum()
    # raw spinning H5 (Mehler with complex index nu = -2i: cos -> cosh 2t; d_mu -> t sinh 2t)
    ch = np.cosh(2*t); sh = np.sinh(2*t)
    H5 = K1*(W*t*sh).sum() + K2*(W*ch).sum() + K3*P0
    return I20, I22, I3, P0, D, H5


def scan(Lam, name, betas, us, BQ, H5mode):
    F, c2, c3 = load_func(name)
    a, c = F['a'], F['c']
    C4 = {+1: 1.0, -1: 41.0}
    best = {'b>=20': (np.inf, None), '5<=b<20': (np.inf, None), 'beta>=2': (np.inf, None), 'beta>=6': (np.inf, None)}
    first_neg = None
    for beta in betas:
        I20, I22, I3, P0, D, H5r = BQ[(Lam, round(beta, 12))]
        H5 = P0 if H5mode == 'pure' else H5r
        bh = Lam*beta
        for h in (+1, -1):
            jmin = 0 if h > 0 else 4
            twoJ1 = Lam*np.sqrt(2*(np.cosh(beta) - 1)/us)
            phys = twoJ1 >= 2*jmin + 1
            B2 = 2*I20 - us*I22
            B3 = h*us*I3
            FW = a*(2*us**2*Lam**2*D - C4[h]*us**3*P0 + h*us**3*H5) - 2*c*us**2*P0
            FF = B2 + B3 + FW
            Nrm = np.abs(B2) + np.abs(B3) + a*(2*us**2*Lam**2*D + 42*us**3*P0 + us**3*abs(H5)) + 2*c*us**2*P0
            rat = np.where(phys, FF/Nrm, np.inf)
            i = int(np.argmin(rat))
            if not np.isfinite(rat[i]):
                continue
            J = 0.5*(twoJ1[i] - 1); m = Lam/np.sqrt(us[i])
            loc = (rat[i], (beta, us[i], h, J, m, bh))
            for k, ok in (('b>=20', bh >= 20), ('5<=b<20', 5 <= bh < 20), ('beta>=2', beta >= 2 and bh >= 20), ('beta>=6', beta >= 6 and bh >= 20)):
                if ok and rat[i] < best[k][0]:
                    best[k] = loc
            if first_neg is None and bh >= 20 and rat[i] < 0:
                first_neg = loc[1] + (rat[i],)
    return best, first_neg


def fmt(b):
    v, (beta, u, h, J, m, bh) = b
    return f"{v:+.4f} (beta {beta:.2f}, {'++' if h > 0 else '+-'}, J {J:.3g}, m/Lam {1/np.sqrt(u):.3g}, bhat {bh:.3g})"


if __name__ == '__main__':
    t0 = time.time()
    pr("# stage5_positivity: leading-Regge model, both forward pieces de-shifted ('pure') vs raw spinning C5 ('raw'); beta <= 30, u in [1e-6, 1]")
    # checks: closed-form cos transform vs GL; Mehler H5 vs mpmath
    F0, c2, c3 = load_func('Q0-full-robust')
    T = [CosT(c2), CosT(np.convolve(c2, [0, 0, 1])), CosT(c3)]
    w = np.array([80., 150., 400., 2000.])
    p, wp = FR.composite_gl(0., 1., 4000, 20)
    ref = np.cos(np.outer(w, p)) @ (wp*np.polynomial.polynomial.polyval(p, c2))
    pr(f"  check: closed-form cos transform vs GL at w = 80, 150, 400, 2000: max abs diff {np.max(np.abs(T[0](w) - ref)):.1e} (values {np.abs(ref).max():.1e})")
    for beta in (0.5, 2.0, 4.0, 8.0):
        x = mp.cosh(beta)
        dP = mp.diff(lambda mu: mp.legenp(mu, 0, x, type=3), mp.mpf(1.5))
        ex = K1*dP + K2*mp.legenp(1.5, 0, x, type=3) + K3*mp.legenp(-0.5, 0, x, type=3)
        t, W = FR.mehler_nodes(beta, 4000)
        H5 = K1*(W*t*np.sinh(2*t)).sum() + K2*(W*np.cosh(2*t)).sum() + K3*W.sum()
        pr(f"  check: raw H5 Mehler vs mpmath at beta = {beta}: {H5:+.10e} vs {float(mp.re(ex)):+.10e} (rel {abs(H5/float(mp.re(ex)) - 1):.1e})")
    pr("  raw spinning H5(beta)/(beta P_{3/2}) and sign change:")
    for beta in (0.5, 2, 4, 6, 7, 8, 12, 20, 30):
        t, W = FR.mehler_nodes(beta, 4000)
        H5 = K1*(W*t*np.sinh(2*t)).sum() + K2*(W*np.cosh(2*t)).sum() + K3*W.sum()
        pr(f"     beta = {beta:5.1f}: H5 = {H5:+.5e}, H5/(beta P_3/2) = {H5/(beta*(W*np.cosh(2*t)).sum()):+.4f}, H5/P0 = {H5/W.sum():+.3e}")
    us = np.concatenate([np.geomspace(1e-6, 0.05, 40), np.linspace(0.06, 1.0, 48)])
    for Lam in (100., 300., 1000.):
        betas = np.unique(np.concatenate([np.geomspace(5/Lam, 1.0, 70), np.linspace(1.0, 12.0, 221), np.linspace(12.0, 30.0, 91)]))
        for name in ('Q0-full-robust', 'Q3-full-delta0.2'):
            ta = time.time()
            Fn, c2n, c3n = load_func(name)
            Tn = [CosT(c2n), CosT(np.convolve(c2n, [0, 0, 1])), CosT(c3n)]
            BQ = {(Lam, round(beta, 12)): beta_quantities(beta, Lam, Tn) for beta in betas}
            for mode in ('pure', 'raw'):
                best, fn = scan(Lam, name, betas, us, BQ, mode)
                pr(f"   Lam = {Lam:5.0f} {name:17s} {mode:4s}: min b>=20: {fmt(best['b>=20'])};  beta>=2: {fmt(best['beta>=2'])};  "
                   f"beta>=6: {fmt(best['beta>=6'])};  [5<=b<20, outside model: {fmt(best['5<=b<20'])}]")
                if fn:
                    beta, u, h, J, m, bh, v = fn
                    pr(f"        first negative (increasing beta): beta {beta:.2f}, {'++' if h > 0 else '+-'}, J {J:.3g}, m/Lam {1/np.sqrt(u):.3g}, value {v:+.3f}")
            pr(f"      [{len(betas)} betas, {time.time() - ta:.0f}s]")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(CACHE, 'stage5_positivity.txt'), 'w').write(buf.getvalue())
