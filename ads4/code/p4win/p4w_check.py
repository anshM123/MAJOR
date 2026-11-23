"""p4win/p4w_check.py (P4_WINDOW.md §3.3): Lemma P4W on the exact scalar.  For the four exact states (p4w_exact) and G = f_ex - f_hyb (the
full beyond-hybrid correction): (a) envelope check |G| <= e_sc (p4w_contour) on every node; (b) the window cost of G against its
series (B2 + S4/m^4 + S6/m^6 + S8/m^8); (c) the complement cost of G (actual, real axis) against the contour bound (uniform and local
rate) with the boundary amplitude taken from the series at nu_w and the envelope e_sc on [p_w, 1]; (d) the same for k = 3 (rho_3:
odd-k interpolation [M] as in sm_exact) with weight u.
-> cache/check.txt, cache/check.json
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
import p4w_common as C, p4w_exact as X, p4w_contour as K

OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))


def main():
    pr(__doc__.split('->')[0])
    res = {}
    for key, r in X.load().items():
        st = r['st']; J = r['J']; m = st['m']; mfl = st['mfl']; beta = st['beta']; u = 1/mfl**2
        n = r['npn']; pn = np.asarray(r['pn']); pw_ = np.asarray(r['pw'])
        nu = np.asarray(r['nus'])[:n]; x = np.asarray(r['x'])[:n]
        E2 = (2 - x**2)/(1 - x**2)**2; E3 = 1/(1 - x**2)**2
        G2 = E2*np.asarray(r['rho'][2])[:n]/2
        G3 = E3*np.asarray(r['rho'][3])[:n]/2 if 3 in r['rho'] else None
        Bud, _ = X.budget(J, mfl)
        env = K.e_sc(beta, mfl)
        ratio = np.max(np.abs(G2)/env(pn))
        P, dPb = C.Pd(beta, nu); P = P.real
        bb = (2*J + 1)/mfl*nu/C.L
        ser = {k: (C.T(k, beta, nu, P, dPb)[0] - nu**(2*k)*C.gflat(k, bb))/m**(2*k) for k in range(1, 5)}
        pr(f"\n=== ({st['Delta']:.0f}, {J}) m/Lam {mfl:.4f} beta {beta:.4f}; budget {Bud:.3e}; envelope check max |G|/e_sc = {ratio:.2f} "
           f"(k=3: {np.max(np.abs(G3)/env(pn)) if G3 is not None else float('nan'):.2f})")
        act2 = lambda f, msk: 2*float((C.psi2(pn)*pw_*f*msk).sum())
        act3 = lambda f, msk: u*float((C.psi3(pn)*pw_*f*msk).sum())
        row = dict(J=J, mfl=mfl, beta=beta, budget=Bud, env_ratio=float(ratio), c={})
        for c in (3, 5, 8):
            nuw = c*m/(2*J + 1); pw0 = nuw/C.L
            win = (nu < nuw).astype(float)
            cw = act2(G2, win); cc = act2(G2, 1 - win)
            sw = [act2(sum(ser[k] for k in range(1, K_ + 1)), win) for K_ in (1, 2, 3, 4)]
            # boundary amplitude from the series (window side), sup over [0.7, 1.3] nu_w, x2 for the travelling-wave split
            nn = np.linspace(0.7*nuw, 1.3*nuw, 61)
            Pn, dPn = C.Pd(beta, nn); bn = (2*J + 1)/mfl*nn/C.L
            Gs = sum((C.T(k, beta, nn, Pn.real, dPn)[0] - nn**(2*k)*C.gflat(k, bn))/m**(2*k) for k in range(1, 3))
            Mw = 2*float(np.max(np.abs(Gs)))
            out = {}
            for mode in ('uniform', 'local'):
                bd, ed, tp = K.contour_bound(beta, mfl, pw0, Mw, env, mode, k=2, Aunif=K.A_curv(beta, J, mfl))
                out[mode] = (bd, ed, tp)
            bd3, ed3, tp3 = K.contour_bound(beta, mfl, pw0, Mw, env, 'local', k=3, u=u, Aunif=K.A_curv(beta, J, mfl))
            cc3 = act3(G3, 1 - win) if G3 is not None else np.nan
            pr(f"  c = {c} (nu_w = {nuw:6.1f}, p_w = {pw0:.4f}): window cost of G {cw/Bud:+.2e} budget; series on window K=1..4: " +
               ' '.join(f'{v/Bud:+.2e}' for v in sw) + f" | complement actual {cc/Bud:+.2e} (k=3: {cc3/Bud:+.2e}) | contour bound "
               f"[bd, edge, top]: uniform [{out['uniform'][0]/Bud:.1e}, {out['uniform'][1]/Bud:.1e}, {out['uniform'][2]/Bud:.0e}], "
               f"local [{out['local'][0]/Bud:.1e}, {out['local'][1]/Bud:.1e}, {out['local'][2]/Bud:.0e}]; k=3 local "
               f"[{bd3/Bud:.1e}, {ed3/Bud:.1e}]  (M_w = {Mw:.2e}, e_sc(1) = {float(env(np.array([1.0]))[0]):.2e})")
            row['c'][c] = dict(cw=cw, cc=cc, cc3=cc3, sw=sw, Mw=Mw, bound={md: list(v) for md, v in out.items()}, bound3=[bd3, ed3, tp3])
        res[str(J)] = row
    json.dump(res, open(os.path.join(C.CACHE, 'check.json'), 'w'), default=float)
    open(os.path.join(C.CACHE, 'check.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
