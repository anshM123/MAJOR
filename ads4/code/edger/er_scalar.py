"""
edger/er_scalar.py (EDGE_R.md §4): the EXACT scalar (Delta_phi = 3, k = 2, Arb Laplace/Barnes via agproof/ag_exact, read-only) at the
gap states of the BINDING Theorem W rows J = 26, 30, 34, 38 (x1 = 0.975-0.964; AG_PROOF had no exact scalar there) and J = 40, 48
(AG_PROOF's), used for three things:
  (1) its (EDGE-R) quantity E_chi[G] = int_{p_w}^1 psi_2 chi G dp, G = ex/2 - E_2 P_J/2 - (P - J0) (P4W units; ag_edge2 convention), at
      the W state's own cutoff and over AG_PROOF's grid, in units S_E2 = 2 V_Fa2(1) e_sc(1);
  (2) its split at the n1 failure point nu_c = 0.94 Delta (x_c = nu_c/m): the part on p <= p_c (where direct graviton data exist) and
      the part on [p_c, 1] (where they do not), and the sign-free size of the latter;
  (3) a test of EXTRAPOLATION protocols that use only data on x <= x_c (as available for the graviton): the edge-window shape model
      Gamma = c(x) D + s(x) x D'/(2J+1) + r(x) (D = P_J(1 - 2x^2), Gamma = m^2 G/(E_2/2); polynomials of degree K in x, fitted on the
      GK15 nodes of [0.8, x1] below x_c plus a 0.0125 grid on [0.8, x_c]) continued to x1, and the amplitude-growth model of HIGHSPIN
      (|Gamma| <= A ((1 - x_l)/(1 - x))^2) used sign-free.  The protocol error on E_chi is measured against (1).
-> cache/scalar.txt, cache/scalar.json
"""
import er_core as R
import os, sys, json, time
import numpy as np
from scipy.special import j0, eval_legendre
import ag_exact as AE           # read-only
import ag_flat as AF            # read-only
import p4w_common as C          # read-only
import p4w_contour as K         # read-only

OUT = []
pr = lambda s='': (print(s, flush=True), OUT.append(s))
GRID_PE = (0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75)
GRID_T = (1.0, 0.8, 0.65, 0.5)
GRID_S = (0.35, 0.45, 0.6)
L = AE.LAM


def carrier(J, x):
    P = eval_legendre(J, 1 - 2*x*x)
    h = 1e-6
    dP = (eval_legendre(J, 1 - 2*(x + h)**2) - eval_legendre(J, 1 - 2*(x - h)**2))/(2*h)
    return P, x*dP/(2*J + 1)


def shape_fit(J, xs, y, K, x0=0.8):
    t = lambda x: (np.asarray(x) - x0)
    D, S = carrier(J, np.asarray(xs))
    A = np.array([t(xs)**k*D for k in range(K + 1)] + [t(xs)**k*S for k in range(K + 1)] + [t(xs)**k for k in range(K + 1)]).T
    sc = np.abs(A).max(0); sc[sc == 0] = 1
    c, *_ = np.linalg.lstsq(A/sc, y, rcond=None); c = c/sc

    def f(x):
        x = np.asarray(x, float); D, S = carrier(J, x)
        B = np.array([t(x)**k*D for k in range(K + 1)] + [t(x)**k*S for k in range(K + 1)] + [t(x)**k for k in range(K + 1)]).T
        return B @ c
    return f


def main(Js):
    t0 = time.time()
    pr(__doc__.split('->')[0])
    res = {}
    for J in Js:
        st = AE.state(J); m = st['m']; beta = st['beta']; mfl = st['mfl']; Dl = st['Delta']
        L3 = AE.build(st, numax=L*1.02, prec=400, conf=(320, 960, 0.5))
        k0 = 0.97*beta*L; T0 = 25.0/k0; pw = 3.0/beta/L
        kmax = (2*J + 1)/mfl/np.sqrt(1 - (1/mfl)**2)
        nseg = int(np.ceil((1 - pw)*1.2*kmax/np.pi*2)) + 40
        xs, ws = np.polynomial.legendre.leggauss(16)
        br = np.linspace(pw, 1.0, nseg + 1)
        p = np.concatenate([a + (b - a)*(xs + 1)/2 for a, b in zip(br[:-1], br[1:])])
        wp = np.concatenate([ws*(b - a)/2 for a, b in zip(br[:-1], br[1:])])
        nu = L*p; x = nu/m
        ex = np.array([float(L3.value_fast(2, float(v)).mid()) for v in nu])/2
        E2h = (2 - x*x)/(1 - x*x)**2/2
        flat = E2h*eval_legendre(J, 1 - 2*x*x)
        cur = C.Pd(beta, nu)[0].real - j0((2*J + 1)/mfl*p)
        G = ex - flat - cur
        ps = K.psi_c(2, p).real
        fl = AF.edge_compare(J, mfl)
        esc1 = float(K.e_sc(beta, mfl)(np.array([1.0]))[0])
        SE2 = 2*fl['VF2']*esc1
        pc = 0.94*Dl/L; xc = pc*L/m; x1 = 1/mfl
        kl = lambda pe: (2*J + 1)/mfl/np.sqrt(1 - (pe/mfl)**2)
        cuts = [(pe, sf*tf*T0) for pe in GRID_PE for tf in GRID_T for sf in GRID_S if kl(pe)*sf*tf*T0 >= 10]
        wr = R.nearest_cut(J, mfl); own = (wr['pe'], wr['sigma'])
        cuts = [own] + cuts
        # exact values at x <= xc on the 'available' node set (0.0125 grid + GK15 of [0.8, x1]) -> shape models
        import hs_core as HC
        xg = np.arange(0.8, xc + 1e-9, 0.0125)
        xe, _, _ = HC.gk15_nodes(0.8, x1)
        xn = np.unique(np.round(np.concatenate([xg, xe[xe <= xc]]), 7))
        Gn = np.array([float(L3.value_fast(2, float(v)).mid()) for v in xn*m])/2
        Pn = eval_legendre(J, 1 - 2*xn*xn); E2n = (2 - xn*xn)/(1 - xn*xn)**2/2
        Gn = Gn - E2n*Pn - (C.Pd(beta, xn*m)[0].real - j0((2*J + 1)*xn))
        Gam_n = m*m*Gn/E2n
        Gam = m*m*G/E2h
        models = {}
        for Kd in (1, 2, 3):
            if len(xn) >= 3*(Kd + 1) + 1:
                f = shape_fit(J, xn, Gam_n, Kd)
                models[f'shape K={Kd}'] = np.where(x > xc, f(x)*E2h/m/m, G)
        sel = x > xc
        last = np.argsort(xn)[-3:]
        A = np.abs(Gam_n[last]).max(); xl = xn.max()
        env = 3*A*np.maximum(1.0, ((1 - xl)/(1 - x))**2)*E2h/m/m     # HIGHSPIN x3 growth envelope (sign-free)
        rows = []
        for (pe, sg) in cuts:
            ch = R.chi(p, pe, sg)
            Ef = float(np.sum(wp*ps*ch*G))
            Ein = float(np.sum((wp*ps*ch*G)[~sel])); Eout = Ef - Ein
            sf_out = float(np.sum((wp*np.abs(ps)*ch*np.abs(G))[sel]))
            cov = float(np.sum((wp*np.abs(ps)*ch*env)[sel]))
            mod = {k: float(np.sum(wp*ps*ch*v)) for k, v in models.items()}
            rows.append(dict(pe=pe, sigma=sg, E=Ef, Ein=Ein, Eout=Eout, sf_out=sf_out, cover=cov, mod=mod))
        r0 = rows[0]
        rr = np.array([abs(r['E']) for r in rows[1:]])/SE2
        pr(f"J = {J} (gap, m/L = {mfl:.4f}, x1 = {x1:.4f}); x_c = {xc:.4f} (nu = 0.94 Delta); S_E2 = {SE2:.3e}; W cutoff (J = {wr['J']}, "
           f"{wr['mfl']:.3f}): p_e = {own[0]}, sigma = {own[1]:.3f}")
        pr(f"   (1) |E_chi[G]| = {abs(r0['E'])/SE2:.3f} S_E2 (own cutoff); AG grid ({len(rows) - 1}): min {rr.min():.3f} median "
           f"{np.median(rr):.3f} max {rr.max():.3f}; max|G| near p = 1: {np.abs(G[p > 0.98]).max():.1f} vs e_sc(1) = {esc1:.1f} "
           f"(rho = max|G|/(E_2/2) = {np.abs(G[sel]/E2h[sel]).max():.3f})")
        pr(f"   (2) split at x_c: E(p <= p_c) = {r0['Ein']/SE2:+.3f}, E(p > p_c) = {r0['Eout']/SE2:+.3f} S_E2; sign-free |.| of the edge part "
           f"{r0['sf_out']/SE2:.2f} S_E2; HIGHSPIN x3 growth envelope (sign-free) {r0['cover']/SE2:.2f} S_E2")
        pr("   (3) extrapolation from x <= x_c (own cutoff; error = model - exact, S_E2): " +
           "; ".join(f"{k}: {(v - r0['E'])/SE2:+.3f}" for k, v in r0['mod'].items()) +
           "  | max over grid: " + "; ".join(f"{k}: {max(abs(r['mod'][k] - r['E']) for r in rows)/SE2:.3f}" for k in r0['mod']))
        res[J] = dict(J=J, mfl=mfl, x1=x1, xc=xc, SE2=SE2, esc1=esc1, rows=rows, nnodes=len(xn), t=time.time() - t0)
        json.dump({str(k): v for k, v in res.items()}, open(os.path.join(R.CACHE, 'scalar.json'), 'w'), default=float)
    pr(f'[{time.time() - t0:.0f}s]')
    open(os.path.join(R.CACHE, 'scalar.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main([int(v) for v in sys.argv[1].split(',')])
