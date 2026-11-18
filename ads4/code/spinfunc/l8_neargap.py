"""
l8_neargap (M6 round 4, doc 22.4): route-A' near-gap positivity with the computed graviton-specific C^_2 correction, per state as the
2x2 (++, +-) matrix (J even >= 4; 1x1 for J = 0, 2), minimum eigenvalue of N^{-1/2} F N^{-1/2}, N = diag(N_s).
  F_hh = base_h (= m2/h2 'base': stripped flat B2 + B3, (P - J0) terms, exact scalar remainder x stripped helicity envelope; no inversion)
         + sum_nodes psi2 w (2-y) Dt_h dc2_h(x; m)/m^2                        (graviton-specific, relative to the Regge point x -> 0)
  F_+- = sum_nodes psi2 w (2-y) offp(x; m)/m                                  (Re part of the O(1/m) off-diagonal; doc 22.3)
  ++ : (2-y) P_J dc2 = (2-y) A+(x),  A+ = P_J d++;
  +- : (2-y) Dt_- dc2 = (2-y)/(1-y)^4 [ A-(x)/m^2 m^2 - m^2 (rho2(x) - d44(x) rho2(0)) ],  A- = d44 Eminus  (rho2 = exact scalar remainder, e1)
m-dependence: model 'asym' = a0(x) + a1(x)/m fitted on Delta = 400...3200 rows; model 'direct' = crystal values interpolated in 1/m through
Delta = 100...3200 (rows with m < 400 re-anchored at x = 0.2 to 'asym', since their x -> 0 extrapolation uses nu <= 20).
x > 0.8: held at the x = 0.8 value (flagged: weight fraction reported).  States with J not in the l8 data (odd J = parity-odd sector,
J > JMAX): left unchanged (listed).  Model 'absmix': as 'asym' but with |off-diagonal| incl. the imaginary part (what doc 21.4 implied).
Also: t* = critical scale factor on the whole correction (both signs) at which n = 1 fails.
Usage: python3 l8_neargap.py [Dg ...] -> cache/l8_neargap_<Dgs>.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0, eval_legendre
from scipy.interpolate import CubicSpline
from spf_common import NearGap, load_states
from common import gl01, Dh_one, SIGMA
from nf_common import flat_parts_nf
from adsaction import envk
import remk
import l7_dc2_analysis as An

XS = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
FIT_D = (400.0, 800.0, 1600.0, 3200.0)
JMAX = int(os.environ.get('L8_JMAX', '10'))      # J = 12: series not usable at m <~ 400 (doc 22.3)


def build_models(store):
    """per J: functions f(x, m, model) -> (A+, A-, O) arrays."""
    mods = {}
    for J, out in store.items():
        if J > JMAX:
            continue
        Ds = sorted(D for D in out if D <= 3200.0)
        if not set(FIT_D) <= set(Ds):
            continue
        PJ = eval_legendre(J, 1 - 2*XS**2); d44 = An.d44(J, 1 - 2*XS**2) if J >= 4 else 0*XS
        rows = {}
        for D in Ds:
            r = out[D]; m2 = r['m']**2
            # exact (not linearised) absolute differences, regular through zeros of P_J, d44 (doc 22.3):
            #   ++: graviton - model = (P_J + rho)(e^dD - 1) -> A+ = P_J (m^2 + s)(exp(d++/m^2) - 1)
            #   +-: graviton (own, rel. to d44 shape) = d44 (1 + s/m^2) e^dD-  -> A- = d44 [(m^2 + s) exp(d+-/m^2) - m^2]
            Ap = PJ*(m2 + r['s'])*np.expm1(r['d++']/m2)
            Am = d44*((m2 + r['s'])*np.exp(r['d+-']/m2) - m2) if J >= 4 else 0*XS
            O = r['offp'] if J >= 4 else 0*XS
            Oa = np.hypot(r['offp'], r['offp_im']) if J >= 4 else 0*XS
            rows[D] = (r['m'], Ap, Am, O, Oa)
        ms = np.array([rows[D][0] for D in FIT_D])
        A = np.array([np.ones_like(ms), 1/ms]).T
        coef = [np.linalg.lstsq(A, np.array([rows[D][k] for D in FIT_D]), rcond=None)[0] for k in (1, 2, 3, 4)]   # each (2, nx)
        resid = [np.abs(A @ c - np.array([rows[D][k] for D in FIT_D])).max() for c, k in zip(coef, (1, 2, 3))]
        def asym(m, coef=coef):
            return [c[0] + c[1]/m for c in coef[:3]]
        def absmix(m, coef=coef):
            return [c[0] + c[1]/m for c in (coef[0], coef[1], coef[3])]
        # direct rows; rows with m < 400 re-anchored at x = 0.2 to 'asym' by a per-state constant (c x LO shape)
        i2 = 2
        drows = []
        for D in Ds:
            m, Ap, Am, O, Oa = rows[D]
            if D < 400.0:
                a = asym(m)
                Ap = Ap + (a[0][i2] - Ap[i2])*PJ/PJ[i2]
                if J >= 4:
                    Am = Am + (a[1][i2] - Am[i2])*d44/d44[i2]
            drows.append((m, [Ap, Am, O]))
        def direct(m, drows=drows):
            im = np.array([1/r[0] for r in drows]); o = np.argsort(im)
            return [np.array([np.interp(1/m, im[o], np.array([drows[i][1][k][ix] for i in o])) for ix in range(len(XS))]) for k in range(3)]
        # grid points too close to zeros of P_J (sign of the scalar remainder lost in the log form) / d44 are not used (doc 22.3)
        okp = np.abs(PJ) >= 0.05; okm = okp & (np.abs(d44) >= 0.02) if J >= 4 else okp
        mods[J] = dict(asym=asym, direct=direct, absmix=absmix, resid=resid, masks=(okp, okm, np.ones_like(okp)))
    return mods


def interp_x(vals, x, masks):
    out = []
    for v, mk in zip(vals, masks):
        xs = np.concatenate([[0.0], XS[mk]])
        cs = CubicSpline(xs, np.concatenate([[0.0], v[mk]]), bc_type='natural')
        out.append(cs(np.minimum(x, XS[mk][-1])))
    return out


def main(Dgs):
    store = pickle.load(open(os.path.join(HERE, 'cache', 'l8_dc2.pkl'), 'rb'))
    mods = build_models(store)
    Jc = sorted(mods)
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    pr(f"# l8_neargap: route-A' near-gap 2x2 with the computed C^_2 graviton-specific correction; computed J = {Jc}")
    pr("#   asym-fit max residuals (c2 units x shape) per J [A+, A-, O]: " +
       "; ".join(f"J={J}: " + "/".join(f"{v:.0f}" for v in mods[J]['resid']) for J in Jc))
    for Dg in Dgs:
        NG = NearGap(Dg); F = NG.F; pn, pw, L = NG.pn, NG.pw, NG.Lam
        blocks = {}
        for h, DF, r0, ref in load_states(Dg):
            mfl, J = r0['mfl'], r0['J']; u = 1/mfl**2; s = SIGMA[h]
            x = pn*np.sqrt(u); y = x*x
            P = remk.Pders_vec(L*pn, r0['eta'], nmax=0)[0]; jz = j0((2*J + 1)*x)
            pf2, wf2 = gl01(max(400, 2*J + 100))
            Js, B2, B3, _ = flat_parts_nf(F, h, u, J, pf2, wf2); i = int(np.searchsorted(Js, J)); B2, B3 = B2[i], B3[i]
            p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
            cv2 = (p2*2*(P - jz)).sum(); cv3 = (p3*u*s*(P - jz)).sum()
            fac = (1 - y)**2 if h == '++' else (1 - y)**-2
            rem = (p2*envk(2, x)*fac*r0['rho'][2]).sum() + (p3*u*s*envk(3, x)*fac*r0['rho'][3]).sum()
            Ns = abs(B2 + cv2) + abs(B3 + cv3)
            blocks.setdefault((round(mfl, 9), J, DF), {})[h] = dict(base=B2 + B3 + cv2 + cv3 + rem, Ns=Ns, r0=r0, x=x, y=y, p2=p2, m=r0['m'])
        for model in ('asym', 'direct', 'absmix'):
            for cover in ('none',):
                recs = []
                wfar = 0.0
                for (mfl, J, DF), b in blocks.items():
                    hs = [h for h in ('++', '+-') if h in b]
                    base = np.array([b[h]['base'] for h in hs]); Ns = np.array([b[h]['Ns'] for h in hs])
                    Jm = J if J in mods else None
                    dF = np.zeros((len(hs), len(hs)))
                    if Jm is not None and (J % 2 == 0 or cover == 'nearest'):
                        e = b[hs[0]]; x, y, p2, m = e['x'], e['y'], e['p2'], e['m']
                        Ap, Am, O = interp_x(mods[Jm][model](m), x, mods[Jm]['masks'])
                        wfar = max(wfar, (p2*(x > XS[-1])).sum()/p2.sum())
                        for k, h in enumerate(hs):
                            if h == '++':
                                dF[k, k] = (p2*(2 - y)*Ap).sum()/m**2
                            elif Jm >= 4:
                                rho = e['r0']['rho'][2] if h in b else 0
                                r0v = np.polyfit(x[:3], b[h]['r0']['rho'][2][:3], 2)[-1]
                                d44 = An.d44(J, 1 - 2*y) if J >= 4 else 1.0
                                dF[k, k] = (p2*(2 - y)/(1 - y)**4*(Am/m**2 - (b[h]['r0']['rho'][2] - d44*r0v))).sum()
                        if len(hs) == 2 and Jm >= 4:
                            dF[0, 1] = dF[1, 0] = (p2*(2 - y)*O).sum()/m
                    recs.append(dict(mfl=mfl, J=J, DF=DF, hs=hs, base=base, Ns=Ns, dF=dF))
                grid = sorted(set(r['mfl'] for r in recs))
                def lam(r, t):
                    Fm = np.diag(r['base']) + t*r['dF']; n = 1/np.sqrt(r['Ns'])
                    return np.linalg.eigvalsh(n[:, None]*Fm*n[None, :])[0]
                def nval(t, delta=0.0):
                    return next((n0 for n0 in grid if all(lam(r, t) > delta for r in recs if r['mfl'] >= n0 - 1e-9)), None)
                w0 = min(recs, key=lambda r: lam(r, 0.0)); w1 = min(recs, key=lambda r: lam(r, 1.0))
                tst = []
                for sgn in (1, -1):
                    lo, hi = 0.0, 1.0
                    while nval(sgn*hi) == grid[0] and hi < 1e4:
                        lo, hi = hi, 2*hi
                    for _ in range(30):
                        mid = 0.5*(lo + hi)
                        lo, hi = (mid, hi) if nval(sgn*mid) == grid[0] else (lo, mid)
                    tst.append(sgn*lo)
                mix = max(recs, key=lambda r: abs(r['dF'][0, 1])/np.sqrt(r['Ns'].prod()) if len(r['hs']) == 2 else 0)
                pr(f"  Delta_g = {Dg} [{model}, uncovered J: {cover}]: n = {nval(1.0)} (no correction: {nval(0.0)}); "
                   f"min lam/N_s = {lam(w1, 1.0):+.4f} [m/Lam={w1['mfl']} J={w1['J']} Df={w1['DF']} {w1['hs']}] "
                   f"(no corr.: {lam(w0, 0.0):+.4f} [m/Lam={w0['mfl']} J={w0['J']}]); t* = ({tst[1]:+.1f}, {tst[0]:+.1f}); "
                   f"max |F_+-|/sqrt(Ns Ns) = {abs(mix['dF'][0, 1])/np.sqrt(mix['Ns'].prod()):.1e} [m/Lam={mix['mfl']} J={mix['J']}]; "
                   f"weight x>0.8: {wfar:.1e}")
                pr("      uncovered (unchanged) states: J = " + str(sorted(set(r['J'] for r in recs if r['J'] not in mods))))
                if model == 'asym' and cover == 'none':
                    dd = [abs(r['dF'][k, k])/r['Ns'][k] for r in recs for k in range(len(r['hs']))]
                    pr(f"      max |dF_hh|/N_s = {max(dd):.2e}, median {np.median(dd):.1e}; worst-state shift of lam_min/N_s: "
                       f"{max(abs(lam(r, 1.0) - lam(r, 0.0)) for r in recs):.2e}")
    sfx = '' if JMAX == 10 else '_Jmax%d' % JMAX
    open(os.path.join(HERE, 'cache', 'l8_neargap_' + '_'.join(map(str, Dgs)) + sfx + '.txt'), 'w').write("\n".join(lines) + "\n")


if __name__ == '__main__':
    main([int(a) for a in sys.argv[1:]] or [100, 300, 1000])
