"""p4win/p4w_exact.py (P4_WINDOW.md §3): the scalar analogue of P4, EXACTLY, on the four Delta_g = 1000 gap states where the exact scalar
(Delta_phi = 3) heavy moments exist on the full nu range (spinfunc/cache/sm_exact_1000.pkl, read-only: Arb Laplace/Barnes route of
exactskel/h13_exact, nodes = Gauss-Legendre in p on [0, 1] + probes).

Per unit leading image (f = moment/2; leading f0 = P):
  f_ex   = ex_2/2                                         exact
  f_hyb  = [E_2(x) P_J(1-2x^2) + 2(P - J0(bhat p))]/2      the scalar hybrid (flat partial wave + leading-Regge curvature), x = nu/m
  B2     = [T_1(nu) - nu^2 g_1(bhat p)]/m^2               beyond-hybrid O(1/m^2) (CHMRS series minus its flat O(y) tower = hybrid's O(y))
  R_ex   = f_ex - f_hyb - B2                              the EXACT P4 remainder (beyond-hybrid, beyond O(1/m^2), all orders)
  R_ser^(K) = sum_{k=2..K} [T_k - nu^{2k} g_k(bhat p)]/m^{2k}    its truncated fixed-beta series (K = 2: P2_CONTRACT's R4)
Costs: k = 2 de-shift action 2 int_0^1 psi2 R dp (Gauss-Legendre nodes of the cache), as fractions of Theorem W's budget H - d2 (bh mode,
p2con/cache/jstar_bh.json) where the state is a Theorem W sample state, else of the hybrid 2x2 budget H recomputed with lsw.floors;
and as K4-equivalents |cost| m^4 / c4a (c4a = |cost(eta^2 nu^4 P)|, P2_CONTRACT's unit).
Window split at nu_w = c m/(2J+1) (c = 1, 3, 5, 10) and nu_w = kappa m (kappa = 0.05, 0.1, 0.2, 0.3): window / complement costs, the
truncation error of the series inside the window, and the complement envelope of R_ex.
-> cache/exact.txt, cache/exact.json
"""
import os, sys, pickle, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0, eval_legendre
import p4w_common as C

OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))
KEYS = [(40, 0.0, 3.0, 'npf2'), (48, 0.0, 3.0, 'd48'), (64, 0.0, 3.0), (100, 0.0, 3.0)]


def load():
    db = pickle.load(open(os.path.join(C.CODE, 'spinfunc', 'cache', 'sm_exact_1000.pkl'), 'rb'))
    return {k: db[k] for k in KEYS}


def budget(J, mfl):
    for s in json.load(open(os.path.join(C.CODE, 'p2con', 'cache', 'jstar_bh.json')))['states']:
        if s['J'] == J and abs(s['mfl'] - mfl) < 1e-6:
            return s['H'] - s['d2'], 'W sample (H - d2, bh)'
    f = C.W.floors(J, mfl)
    fl = [v['floor'] for v in f.values()]
    H = fl[0]*fl[1]/(fl[0] + fl[1]) if len(fl) == 2 else fl[0]
    return H, 'H (recomputed, d2 not available)'


def main():
    recs = load()
    res = {}
    pr(__doc__.split('->')[0])
    for key, r in recs.items():
        st = r['st']; J = r['J']; m = st['m']; mfl = st['mfl']; beta = st['beta']
        n = r['npn']; pn = np.asarray(r['pn']); pw = np.asarray(r['pw'])
        nu_all = np.asarray(r['nus']); x_all = np.asarray(r['x'])
        bh = (2*J + 1)/mfl
        P, dPb = C.Pd(beta, nu_all)
        P = P.real
        Tk = {}; cond = {}
        for k in range(1, 5):
            Tk[k], cond[k] = C.T(k, beta, nu_all, P, dPb)
        bb = bh*nu_all/C.L
        g = {k: C.gflat(k, bb) for k in range(1, 5)}
        ex = np.asarray(r['ex'][2])
        E2 = (2 - x_all**2)/(1 - x_all**2)**2
        f_ex = ex/2
        f_hyb = (E2*eval_legendre(J, 1 - 2*x_all**2) + 2*(P - j0(bb)))/2
        rho_chk = np.max(np.abs((ex - 2*f_hyb)/E2 - np.asarray(r['rho'][2])))
        B2 = (Tk[1] - nu_all**2*g[1])/m**2
        R_ex = f_ex - f_hyb - B2
        Rk = {k: (Tk[k] - nu_all**(2*k)*g[k])/m**(2*k) for k in range(2, 5)}
        Rser = {K: sum(Rk[k] for k in range(2, K + 1)) for K in (2, 3, 4)}
        # truncation of the full (not hybrid-subtracted) series: f_ex - (P + sum_k T_k/m^2k)
        S = {0: P}
        for K in range(1, 5):
            S[K] = S[K - 1] + Tk[K]/m**(2*K)
        psw = C.psi2(pn)*pw
        act = lambda f: 2*float((psw*f[:n]).sum())
        c4a = abs(act(np.cosh(beta)**2*nu_all**4*P))
        Bud, budname = budget(J, mfl)
        nsum_err = max(abs(2*float((psw*(np.asarray(ra[2])[:n] - np.asarray(r['rho'][2])[:n])*E2[:n]/2).sum())) for ra in r['rho_alts']) if r['rho_alts'] else 0.0
        pr(f"\n=== exact scalar, (Delta, J) = ({st['Delta']:.0f}, {J}), m/Lam = {mfl:.4f}, beta = {beta:.4f}, bhat = {bh:.1f}; budget {Bud:.3e} [{budname}]")
        pr(f"   checks: rho recomputed vs cache {rho_chk:.1e}; n-sum spread of the cost of R_ex: {nsum_err:.1e} (= {nsum_err/Bud:.1e} budget)")
        row = dict(J=J, mfl=mfl, m=m, beta=beta, budget=Bud, c4a=c4a)
        # ---- full-range costs
        cR = act(R_ex); cB2 = act(B2)
        pr(f"   full range: cost(R_ex) = {cR:+.3e} = {cR/Bud:+.2e} budget (K4-equiv {abs(cR)*m**4/c4a:.3g}); cost(B2) = {cB2/Bud:+.2e} budget")
        for K in (2, 3, 4):
            c = act(Rser[K])
            pr(f"   full range: truncated series R_ser^({K}) cost = {c/Bud:+.2e} budget (K4-equiv {abs(c)*m**4/c4a:.3g})")
            row[f'Kser{K}'] = abs(c)*m**4/c4a
        row.update(cost_ex=cR, K_ex=abs(cR)*m**4/c4a, cost_B2=cB2)
        # ---- truncation of the series vs nu (pointwise, relative to the leading P envelope ~ 1)
        pr("   pointwise |f_ex - S_K| (K = 1..4) and |R_ex - R_ser^(K)| (K = 2..4) at selected nu (GL nodes + probes):")
        idx = np.argsort(nu_all)
        sel = [i for i in idx if nu_all[i] < 0.95*C.L]
        picks = [sel[int(t*(len(sel) - 1))] for t in np.linspace(0, 1, 18)]
        for i in sorted(set(picks), key=lambda i: nu_all[i]):
            pr(f"     nu = {nu_all[i]:7.2f} (nu beta = {nu_all[i]*beta:6.2f}, x = {x_all[i]:.3f}): |f_ex-S_K| = " +
               ' '.join(f'{abs(f_ex[i] - S[K][i]):.1e}' for K in range(1, 5)) +
               f" | R_ex = {R_ex[i]:+.2e}, |R_ex - R_ser^K| = " + ' '.join(f'{abs(R_ex[i] - Rser[K][i]):.1e}' for K in (2, 3, 4)) +
               f"  (cond T4 {cond[4][i]/max(abs(Tk[4][i]), 1e-300):.1e})")
        # ---- window splits
        pr("   window splits: cost of R_ex on nu < nu_w (window) / nu >= nu_w (complement); series R_ser^(2,4) on the window; "
           "max|R_ex - R_ser^(4)|/max|R_ser^(4)| on the window")
        wins = [('c', c, c*m/(2*J + 1)) for c in (1, 3, 5, 10)] + [('kappa', kp, kp*m) for kp in (0.05, 0.1, 0.2, 0.3)]
        row['win'] = []
        for kind, par, nuw in wins:
            msk = (nu_all[:n] < nuw).astype(float)
            mfull = np.concatenate([msk, np.zeros(len(nu_all) - n)])
            cw = act(R_ex*mfull); cc = act(R_ex*(1 - mfull))
            cs2 = act(Rser[2]*mfull); cs4 = act(Rser[4]*mfull)
            inw = nu_all < nuw
            tr = (np.max(np.abs(R_ex[inw] - Rser[4][inw]))/max(np.max(np.abs(Rser[4][inw])), 1e-300)) if inw.any() else np.nan
            tr2 = (np.max(np.abs(R_ex[inw] - Rser[2][inw]))/max(np.max(np.abs(Rser[2][inw])), 1e-300)) if inw.any() else np.nan
            env = np.max(np.abs(R_ex[~inw & (nu_all < 0.99*C.L)])) if (~inw).any() else np.nan
            pr(f"     {kind} = {par:5g} (nu_w = {nuw:7.2f}, {int(msk.sum())} nodes): window {cw/Bud:+.2e}, complement {cc/Bud:+.2e} budget | "
               f"series on window: R4 {cs2/Bud:+.2e}, R4+R6+R8 {cs4/Bud:+.2e} | rel. trunc. err on window K=2: {tr2:.1e}, K=4: {tr:.1e} | "
               f"max|R_ex| on complement {env:.2e}")
            row['win'].append(dict(kind=kind, par=par, nuw=nuw, cw=cw, cc=cc, cs2=cs2, cs4=cs4, tr2=tr2, tr4=tr, env=env))
        # ---- complement envelope vs x (for the AG amplitude premise): max |R_ex| in x-bins, in units of m^-4 nu^4-ish and absolute
        pr("   complement amplitude of R_ex (sup over x-bins; GL nodes): x-bin: max|R_ex|, max|R_ex| (1-x^2)^3, max|R_ex| m^4/nu^4")
        xs = x_all[:n]; Rn = R_ex[:n]; nun = nu_all[:n]
        row['envbins'] = []
        for a, b in ((0.0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 0.9), (0.9, 0.97)):
            s_ = (xs >= a) & (xs < b)
            if s_.any():
                M0 = np.max(np.abs(Rn[s_])); M3 = np.max(np.abs(Rn[s_])*(1 - xs[s_]**2)**3); M4 = np.max(np.abs(Rn[s_])*m**4/nun[s_]**4)
                pr(f"     x in [{a:.2f}, {b:.2f}): {M0:.2e}  {M3:.2e}  {M4:.2e}")
                row['envbins'].append((a, b, M0, M3, M4))
        res['%d' % J] = row
    json.dump(res, open(os.path.join(C.CACHE, 'exact.json'), 'w'), default=float)
    open(os.path.join(C.CACHE, 'exact.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
