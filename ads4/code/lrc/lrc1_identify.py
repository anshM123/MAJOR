r"""
lrc1_identify (LRC_CLOSE.md §1): what the exact (all-orders-in-rho) CLS-D per-nu kernel IS.

(a) The O(rho) deviation of the nu-dependent bases, as LEMMAS_CLOSE §3.3 measured it: re-run lemmas/lem3_k2nu_lrc.part_b (imported,
    read-only) so the starting point is reproduced.
(b) The exact scalar j'=0 per-nu kernel in the H7 §18.1 Regge frame (x1 = -y, x3 = y, x2 = -x4 = x/x^2, x0 = |y| n_theta on the circle cut
    out by the two light-ray delta functions, z0 on the celestial circle), with the nu-dependence and the (x0, z0)-dependent nu-independent
    factors read off ads4/code/repro_engine/integrand.py (header; CLS-D eq. spinningfunctional_definition):
        K_nu(r, beta) = int dtheta int dphi  A^{i nu/2} |Nn|^{-1/2} m24^{3/4} (m20 m40)^{-1/4} c13^{5/4} (c10 c30)^{-1/4},
        A = Nn^2 c13/(m24 m40 m20 c10 c30),  Nn = 2(z0.x40 (-m20) - z0.x20 (-m40)), c13 = 2 x10.x30, c10 = 2 z0.x10, c30 = -2 z0.x30.
    The dual factor equals V_{0,42}^{l} (X24/(X20 X40))^{(1+l)/2}, l = -1/2 + i nu: the 3-point structure of a (Delta, J) = (1, l) operator with
    weight-0 partners; the celestial factor its (13) partner.  So K_nu should be a (13)(42)-channel Casimir eigenfunction with (Dp, J') = (1, l)
    (or its shadow Dp = 2), times a nu-INDEPENDENT prefactor.
    Test (prefactor-free): Q(x) = [K_nu(x)/K_nu0(x)] / [Gh_nu(x)/Gh_nu0(x)] must be independent of x = (rho, beta) for the right Dp.
Output: cache/lrc1_identify.txt.  ~ 1-2 min.
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'lemmas'))
import numpy as np
import lrc_blocks as L

OUT = []
T0 = time.time()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def mdot(u, v):
    return -u[..., 0]*v[..., 0] + u[..., 1]*v[..., 1] + u[..., 2]*v[..., 2]


def grid_u(ds=0.02, smax=60.0):
    """logistic map u = pi/(1+e^{-s}) on (0, pi): u ~ pi e^s at 0, pi - u ~ pi e^{-s} at pi; returns u, du/ds * ds."""
    s = np.arange(-smax, smax + ds/2, ds)
    u = np.pi/(1 + np.exp(-s))
    w = u*(1 - u/np.pi)*ds
    return u, w


def kernel_parts(r, beta, nth=128, ds=0.02):
    """log of the nu-independent factor (real, all factors positive) and log A on the (theta, u) grid, plus quadrature weights."""
    th = 2*np.pi*np.arange(nth)/nth
    u1, w1 = grid_u(ds)
    u = np.concatenate([u1, u1 + np.pi]); wu = np.concatenate([w1, w1])
    TH = th[:, None]; PH = TH + u[None, :]
    xh = np.array([np.cosh(beta), np.sinh(beta), 0.0]); yh = np.array([1.0, 0.0, 0.0])
    x2 = -xh; x4 = xh; x1 = -r*yh; x3 = r*yh                # |x| = 1, |y| = r
    n = np.stack([0*TH + 0*PH, np.cos(TH) + 0*PH, np.sin(TH) + 0*PH], -1)
    x0 = r*n
    z0 = np.stack([np.ones_like(PH), np.cos(PH), np.sin(PH)], -1)
    x10 = x1 - x0; x30 = x3 - x0; x20 = x2 - x0; x40 = x4 - x0
    m20 = -mdot(x20, x20); m40 = -mdot(x40, x40); m24 = -mdot(x2 - x4, x2 - x4)
    Nn = 2*(mdot(z0, x40)*(-m20) - mdot(z0, x20)*(-m40))
    c13 = 2*mdot(x10, x30)
    U = 0*TH + u[None, :]
    c10 = 4*r*np.sin(U/2)**2; c30 = 4*r*np.cos(U/2)**2          # = 2 z0.x10, -2 z0.x30 exactly (no cancellation at u -> 0, pi)
    A = Nn**2*c13/(m24*m40*m20*c10*c30)
    signs = (np.all(m20 > 0), np.all(m40 > 0), m24 > 0, np.all(c13 > 0), np.all(c10 > 0), np.all(c30 > 0), np.all(Nn > 0) or np.all(Nn < 0), np.all(A > 0))
    logF = (-0.5*np.log(np.abs(Nn)) + 0.75*np.log(m24) - 0.25*np.log(m20*m40) + 1.25*np.log(c13) - 0.25*np.log(c10*c30))
    W = (2*np.pi/nth)*wu[None, :]
    return logF, np.log(A), W, signs


def K_exact(parts, nus):
    logF, logA, W, _ = parts
    out = []
    for nu in nus:
        out.append(np.sum(W*np.exp(logF + 0.5j*nu*logA)))
    return np.array(out)


def Gh_at(Dp, nu, rho, beta, N=160):
    B, m = L.block_coeffs(Dp, nu, N)
    Pt = L.legendre_scaled(nu, beta, N + 1)
    Y = L.radial_Y(B, m, Pt, beta)
    return L.Gh(Y, rho)[0]


def part_a():
    pr('(a) LEMMAS_CLOSE §3.3 (b) reproduced: relative x0-corrections of the nu-dependent bases are rho = r e^B, not r')
    import lem3_k2nu_lrc as LEM
    LEM.pr = pr                     # route its printing into this log
    LEM.part_b()


def part_b():
    pr('(b) exact scalar j\'=0 per-nu kernel (H7 frame, 2-d quadrature over x0-circle and z0) vs continuous-spin blocks Gh_{Dp, -1/2+i nu}')
    nus = np.array([0.4, 0.8, 1.3, 2.7])
    nu0 = 0.8
    pts = [(rho, beta) for rho in (0.25, 0.55, 0.85) for beta in (0.6, 1.4, 2.2)]
    # quadrature convergence check at one point
    p1 = kernel_parts(0.55*np.exp(-2.2), 2.2, nth=128, ds=0.02)
    p2 = kernel_parts(0.55*np.exp(-2.2), 2.2, nth=192, ds=0.012)
    k1 = K_exact(p1, nus); k2 = K_exact(p2, nus)
    pr(f'    quadrature check (rho = 0.55, beta = 2.2): max rel diff between grids = {np.max(np.abs(k1/k2 - 1)):.1e};  '
       f'all bases of fixed sign, A > 0: {p1[3]}')
    rows = []
    for (rho, beta) in pts:
        r = rho*np.exp(-beta)
        K = K_exact(kernel_parts(r, beta), nus)
        Kr = K/K[list(nus).index(nu0)]
        G = {Dp: np.array([Gh_at(Dp, nu, rho, beta) for nu in nus]) for Dp in (1.0, 2.0)}
        Q = {Dp: Kr/(G[Dp]/G[Dp][list(nus).index(nu0)]) for Dp in (1.0, 2.0)}
        rows.append((rho, beta, Kr, Q, r, K[list(nus).index(nu0)], G[1.0][list(nus).index(nu0)]))
        pr(f'    rho = {rho:.2f}, beta = {beta:.1f} (r = {r:.3e}):  '
           + '  '.join(f'Q_{{Dp=1}}(nu={nu}) = {Q[1.0][i].real:+.8f}{Q[1.0][i].imag:+.1e}i' for i, nu in enumerate(nus) if nu != nu0))
    for Dp in (1.0, 2.0):
        for i, nu in enumerate(nus):
            if nu == nu0:
                continue
            qs = np.array([row[3][Dp][i] for row in rows])
            pr(f'    Dp = {Dp:.0f}, nu = {nu}: Q over the 9 points: mean {np.mean(qs.real):+.8f}, max rel spread {np.max(np.abs(qs/np.mean(qs) - 1)):.2e}')
    # the leading-Regge-only comparison (Gh -> P_l): how much nu-dependent structure beyond LR the kernel carries
    for i, nu in enumerate(nus):
        if nu == nu0:
            continue
        qs = []
        for (rho, beta, Kr, Q, r_, k0, g0) in rows:
            import mpmath as mp
            pl = complex(mp.legenp(-0.5 + 1j*nu, 0, mp.cosh(beta), type=3)); pl0 = complex(mp.legenp(-0.5 + 1j*nu0, 0, mp.cosh(beta), type=3))
            qs.append(Kr[i]/(pl/pl0))
        qs = np.array(qs)
        pr(f'    control (LR only, Gh -> P_l): nu = {nu}: max rel spread {np.max(np.abs(qs/np.mean(qs) - 1)):.2e}')
    # the nu-independent prefactor pref(x) = K_nu0/Gh_nu0: a pure power of r (standard 4-pt prefactor with equal weights) ?
    lr = np.log([row[4] for row in rows]); pf = np.array([row[5]/row[6] for row in rows])
    ph = np.angle(pf); lp = np.log(np.abs(pf))
    A = np.vstack([lr, np.ones_like(lr)]).T
    coef, *_ = np.linalg.lstsq(A, lp, rcond=None)
    pr(f'    prefactor K_nu0/Gh_nu0 over the 9 points: phase spread {np.ptp(ph):.1e} rad; log|pref| = {coef[0]:.6f} log r + {coef[1]:.4f}, '
       f'max residual {np.max(np.abs(A @ coef - lp)):.1e}  (a pure power of r: the exact f-stripped kernel is r^p x Gh_(1, l), no other x-dependence)')


if __name__ == '__main__':
    pr('# lrc1_identify: the exact CLS-D per-nu kernel is a continuous-spin (13)(42) Casimir eigenfunction (LRC_CLOSE.md §1)')
    part_a(); pr(f'  [{time.time() - T0:.0f}s]')
    part_b(); pr(f'  [{time.time() - T0:.0f}s]')
    with open(os.path.join(HERE, 'cache', 'lrc1_identify.txt'), 'w') as f:
        f.write('\n'.join(OUT) + '\n')
