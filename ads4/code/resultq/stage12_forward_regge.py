"""Result Q, Stage 1 (RESULT_Q_STAGE12.md §1.4): does the flat '2JJ growth' of -d_{p^2}B4|_0 survive AdS curvature?

Leading-Regge, pure-moment (de-shifted) model of the AdS lift: a pure spin-k Regge moment normalised to the flat B_k heavy action acts on
a heavy block (m, J) as  (flat normal form) x P_{-1/2+i nu}(cosh beta),  cosh beta = 1 + (2J+1)^2/(2 m^2)  (CHMRS; EE: same in both
helicity channels, no off-diagonal, EE_FUNCTIONAL.md §4).  Forward pieces at nu = 0 (flat p = nu/Lam):
   B4|_0      -> 2u^2 P0(beta),                      P0 := P_{-1/2}(cosh beta)
   dB4|_0     -> Lam^2 d_{nu^2}[2u^2 P_{-1/2+i nu}]_0 + c4_h u^3 P0 = -2 u^2 Lam^2 D(beta) + c4_h u^3 P0   (c4_h: flat constants, kept as model)
   B5|_0      -> sigma_h u^3 P0
with D(beta) := -d_{nu^2} P_{-1/2+i nu}(cosh beta)|_{nu=0}.  Mehler-Dirichlet (standard):
   P_{-1/2+i nu}(cosh beta) = (sqrt 2/pi) int_0^beta cos(nu t) (cosh beta - cosh t)^{-1/2} dt
=> D(beta) = (sqrt 2/(2 pi)) int_0^beta t^2 (cosh beta - cosh t)^{-1/2} dt > 0 for every beta > 0   [P]
   R(beta) := D/P0 = <t^2>_w / 2 with w = (cosh beta - cosh t)^{-1/2} on [0, beta]; numerically R/(beta^2/4) = 1 - O(beta^2), 0.91 (beta=4), 0.75 (beta=20).
Flat limit: Lam^2 R -> Lam^2 beta^2/4 = u (J + 1/2)^2 ~ JJ u  (the flat 'JJ u' growth).

Part (b): the full Stage-0 functionals (cache/functionals_q.json) in this leading-Regge model on beta in [beta_min, 12], u in (0, 1],
both channels, Lam = 100, 300, 1000:
   FF_R = int psi2 (2 - u p^2) P_{-1/2+i Lam p} dp + sigma u int psi3 P_{-1/2+i Lam p} dp
          + a [2 u^2 Lam^2 D - c4_h u^3 P0 + sigma u^3 P0] - 2 c u^2 P0           (c = 1/2)
This is a MODEL ([M]): flat finite-b corrections and beyond-leading-Regge corrections are dropped; it tests only the large-beta (Regge)
mechanism N3 of RESULT_Q_STAGE0.md §2 in AdS.  Output: cache/stage12_forward_regge.txt.  ~1-2 min, 2 cores.
"""
import os, io, json, time
os.environ.setdefault('OMP_NUM_THREADS', '2')
import numpy as np
import mpmath as mp
from functools import lru_cache


@lru_cache(maxsize=None)
def _gl(n):
    return np.polynomial.legendre.leggauss(n)


def composite_gl(a, b, panels, order=20):
    x, w = _gl(order)
    edges = np.linspace(a, b, panels + 1)
    h = np.diff(edges)
    pts = (edges[:-1, None] + 0.5 * h[:, None] * (x[None, :] + 1)).ravel()
    wts = (0.5 * h[:, None] * w[None, :]).ravel()
    return pts, wts

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def mehler_nodes(beta, n=1200):
    """nodes t_i and weights W_i with sum W_i g(t_i) = (sqrt2/pi) int_0^beta g(t) (cosh beta - cosh t)^{-1/2} dt
    (substitution t = beta (1 - s^2), removes the endpoint singularity)."""
    s, ws = composite_gl(0., 1., max(1, n // 20), 20)
    t = beta * (1 - s ** 2)
    dt = 2 * beta * s
    den = np.sqrt(np.cosh(beta) - np.cosh(t))
    # cosh beta - cosh t = 2 sinh((beta+t)/2) sinh((beta-t)/2), stable form
    den = np.sqrt(2 * np.sinh((beta + t) / 2) * np.sinh((beta - t) / 2))
    W = np.sqrt(2) / np.pi * ws * dt / den
    return t, W


def P0_D(beta, n=1200):
    t, W = mehler_nodes(beta, n)
    return W.sum(), 0.5 * (W * t ** 2).sum()


def check_mehler():
    worst = 0.
    for beta in (0.01, 0.3, 1.0, 3.0, 8.0):
        for nu in (0.0, 0.7, 5.0):
            t, W = mehler_nodes(beta)
            v = (W * np.cos(nu * t)).sum()
            ref = float(mp.re(mp.legenp(-0.5 + 1j * nu, 0, mp.cosh(beta), type=3)))
            worst = max(worst, abs(v / ref - 1))
            # derivative check for D at nu = 0 (central difference in nu^2)
    h = 1e-3
    dd = []
    for beta in (0.3, 2.0, 6.0):
        P0, D = P0_D(beta)
        fd = -(float(mp.re(mp.legenp(-0.5 + 1j * h, 0, mp.cosh(beta), type=3))) - float(mp.re(mp.legenp(-0.5, 0, mp.cosh(beta), type=3)))) / h ** 2
        dd.append(abs(D / fd - 1))
    return worst, max(dd)


def load_psi(name):
    F = json.load(open(os.path.join(HERE, 'cache', 'functionals_q.json')))[name]
    return F, np.array(F['psi2'], float), np.array(F['psi3'], float) if len(F['psi3']) else np.zeros(1)


def cos_transforms(c2, c3, omegas):
    """Psi2hat0(w) = int psi2 cos(w p), Psi2hat2 = int psi2 p^2 cos(w p), Psi3hat = int psi3 cos(w p), p in [0,1]
    (composite Gauss-Legendre, 20 points per panel, panel width <= 4/omega_max)."""
    p, wp = composite_gl(0., 1., int(max(omegas.max(), 50.) / 4) + 10, 20)
    ps2 = np.polynomial.polynomial.polyval(p, c2); ps3 = np.polynomial.polynomial.polyval(p, c3)
    out = np.zeros((3, len(omegas)))
    for i0 in range(0, len(omegas), 512):
        om = omegas[i0:i0 + 512]
        C = np.cos(np.outer(om, p))
        out[0, i0:i0 + 512] = C @ (wp * ps2)
        out[1, i0:i0 + 512] = C @ (wp * ps2 * p ** 2)
        out[2, i0:i0 + 512] = C @ (wp * ps3)
    return out


def main():
    t0 = time.time()
    e1, e2 = check_mehler()
    pr(f"# stage12_forward_regge.  Mehler-Dirichlet quadrature vs mpmath legenp: max rel err {e1:.1e}; D vs finite difference: {e2:.1e}")
    pr("(a) D(beta) = -d_{nu^2}P_{-1/2+i nu}(cosh beta)|_0 > 0  [P: positive integrand];  R = D/P0;  flat value beta^2/4 ~ (2J+1)^2/(4m^2)")
    pr("    beta      P0           D            R          R/(beta^2/4)   R/(beta^2/2)")
    for beta in (1e-3, 0.01, 0.1, 0.3, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0):
        P0, D = P0_D(beta, 2000)
        R = D / P0
        pr(f"   {beta:7.3f}  {P0:.6e}  {D:.6e}  {R:.6e}  {R / (beta ** 2 / 4):.6f}      {R / (beta ** 2 / 2):.6f}")
    pr("    => the forward d_{p^2} action -2u^2 Lam^2 D(beta) has the flat sign (growth) for ALL beta, in both helicity channels; relative to")
    pr("       B4|_0 it grows like Lam^2 R(beta): = JJ u in the flat region, then (0.75-1) x Lam^2 beta^2/4 (log^2 J growth) for beta >~ 1 (table). [P] (model: pure moment)")
    # (b) full functionals
    C4 = {+1: 1.0, -1: 41.0}
    names = ['Q0-full-robust', 'Q0-min-robust (no B3)', 'Q0-r5-psi3zero-robust (AdS-style weights)']
    betas = np.concatenate([np.geomspace(2e-3, 0.05, 25), np.linspace(0.06, 12, 140)])
    us = np.concatenate([np.geomspace(1e-4, 0.05, 15), np.linspace(0.06, 1.0, 40)])
    pr("(b) leading-Regge pure-moment model of the Stage-0 functionals (c = 1/2), min over u in [1e-4, 1], beta in [2e-3, 12], both channels.")
    pr("    physical states only (b = Lam beta >= (2 J_min + 1) sqrt(u), J_min = 0 (++), 4 (+-)); minima of FF_R/Nrm over: all, b >= 20, beta >= 0.1, beta >= 1")
    pr("    Nrm = |B2 part| + |B3 part| + a*(2u^2 Lam^2 D + 42 u^3 P0) + 2c u^2 P0")
    for Lam in (100., 300., 1000.):
        for nm in names:
            F, c2, c3 = load_psi(nm)
            a, c = F['a'], F['c']
            omax = Lam * betas.max() + 10
            om = np.linspace(0, omax, int(omax * 8) + 2)
            T = cos_transforms(c2, c3, om)
            thr = {'phys': lambda b: b > 0, 'b>=20': lambda b: Lam * b >= 20, 'beta>=0.1': lambda b: b >= 0.1, 'beta>=1': lambda b: b >= 1.0}
            best = {k: (np.inf, None) for k in thr}
            for beta in betas:
                t, W = mehler_nodes(beta, 1600 if Lam * beta > 2000 else 800)
                tr = [np.interp(Lam * t, om, T[i]) for i in range(3)]
                I20 = (W * tr[0]).sum(); I22 = (W * tr[1]).sum(); I3 = (W * tr[2]).sum()
                P0 = W.sum(); D = 0.5 * (W * t ** 2).sum()
                for h in (+1, -1):
                    jmin = 0 if h > 0 else 4
                    phys = Lam * beta >= (2 * jmin + 1) * np.sqrt(us)        # b = (2J+1) sqrt(u) >= (2 J_min + 1) sqrt(u)
                    B2 = 2 * I20 - us * I22
                    B3 = h * us * I3
                    FW = a * (2 * us ** 2 * Lam ** 2 * D - C4[h] * us ** 3 * P0 + h * us ** 3 * P0) - 2 * c * us ** 2 * P0
                    FF = B2 + B3 + FW
                    Nrm = np.abs(B2) + np.abs(B3) + a * (2 * us ** 2 * Lam ** 2 * D + 42 * us ** 3 * P0) + 2 * c * us ** 2 * P0
                    rat = np.where(phys, FF / Nrm, np.inf)
                    i = int(np.argmin(rat))
                    for k, fn in thr.items():
                        if fn(beta) and rat[i] < best[k][0]:
                            best[k] = (rat[i], (beta, us[i], h))
            pr(f"    Lam = {Lam:6.0f}  {nm:45s} " + ";  ".join(
                f"{k}: {best[k][0]:+.4f} (beta={best[k][1][0]:.3f}, u={best[k][1][1]:.4f}, {'++' if best[k][1][2] > 0 else '+-'})" for k in thr))
    pr(f"(done in {time.time() - t0:.0f}s)")
    open(os.path.join(CACHE, 'stage12_forward_regge.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
