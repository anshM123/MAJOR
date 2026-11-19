"""
n2_neargap_split (doc 24.8): copy of l9_neargap with the graviton correction on nodes x < N2_XC scaled by TLOW[0] (tolerance of the
uncontrolled small-x region).  Generated from l9_neargap.py by sed (see run_all.sh); do not edit by hand.

l9_neargap (M6 round 5, doc 23.4): route-A' near-gap positivity per state as the 2x2 (++, +-) matrix with the computed graviton-specific
corrections of C^_2 and (optionally) C^_3, from the LINEAR-form grids of l9_analysis (nu-only factor removed order by order), with DIRECT
values at every computed Delta row (no extrapolation in m inside the grid; linear in 1/m between rows; rows: K3 where computed, else K2).
  F_hh = base_h (as l8_neargap/h2: stripped flat B2 + B3, (P - J0) terms, exact scalar remainder x envelopes; no inversion)
       + sum_nodes psi2 w env_{2,h} [Gamma_{2,h}(x; m)/m^2 - (rho2 - D_h rho2(0))]
       + sum_nodes psi3 w u env_{3,h} [Gamma_{3,h}(x; m)/m^2 - sigma_h (rho3 - D_h rho2(0))]          (if C^_3 included)
  F_+- = sum_nodes psi2 w env_{2,off} Gamma_{2,off}/m (+ psi3 u env_{3,off} Gamma_{3,off}/m)          (Re parts only)
  (graviton normalised to its k = 2 Regge value per helicity, model to its own k = 2 Regge value: the common per-state constant is a
   congruence when both k are corrected, so n does not depend on it; with C^_2 alone it is the doc-22 convention.)
Rows above the largest computed m: held at the largest row.  x beyond a row's last x: held.  J not covered: unchanged (listed).
t*: scale factor on the whole correction at which n = 1 first fails.
Usage: python3 l9_neargap.py variant Dg [Dg ...]   variant in {C2_K2, C2_K3, C23H_K2, C23H_K3, C23V_K3}
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
from scipy.interpolate import CubicSpline
from spf_common import NearGap, load_states
from common import gl01, SIGMA
from nf_common import flat_parts_nf
from adsaction import envk
import remk
import l9_analysis as A

XC = float(os.environ.get("N2_XC", "0.2")); TLOW = [1.0]
JMAX = int(os.environ.get('L9_JMAX', '10'))


def rows_for(J, K, kinds):
    """merged rows: K3 rows (if K == 3) override K2 rows at the same Delta."""
    g = A.grid(J, 2, kinds)
    if K == 3:
        try:
            g3 = A.grid(J, 3, kinds)
        except Exception:
            g3 = {}
        for D, r in g3.items():
            r['K'] = 3; g[D] = r
    return {D: r for D, r in g.items() if D <= 3200.0}


def model(rows, J, kinds):
    Ds = sorted(rows); ms = np.array([rows[D]['m'] for D in Ds])
    keys = [(kd, e) for kd in kinds for e in (('++', '+-', 'off') if J >= 4 else ('++',))]
    spl = {D: {k: CubicSpline(rows[D]['xs'], rows[D][k], bc_type='natural') for k in keys} for D in Ds}
    xmax = {D: rows[D]['xs'][-1] for D in Ds}
    def f(x, m):
        im = 1/m; w = {}
        if m <= ms[0]:
            i0 = i1 = 0; t = 0.0
        elif m >= ms[-1]:
            i0 = i1 = len(Ds) - 1; t = 0.0
        else:
            i1 = int(np.searchsorted(ms, m)); i0 = i1 - 1
            t = (1/ms[i0] - im)/(1/ms[i0] - 1/ms[i1])
        out = {}
        for k in keys:
            v0 = spl[Ds[i0]][k](np.minimum(x, xmax[Ds[i0]])); v1 = spl[Ds[i1]][k](np.minimum(x, xmax[Ds[i1]]))
            out[k] = (1 - t)*v0 + t*v1
        return out
    return f, Ds


def main(variant, Dgs):
    kinds = ('C2',) if variant.startswith('C2_') else ('C2', 'C3H' if 'H' in variant.split('_')[0] else 'C3V')
    K = int(variant[-1])
    mods = {}
    for J in range(0, JMAX + 1, 2):
        try:
            rows = rows_for(J, K, kinds)
        except Exception as ex:
            print(f"J={J}: no grid ({ex})"); continue
        if len(rows) >= 2:
            mods[J] = model(rows, J, kinds) + ({D: rows[D].get('K', 2) for D in rows},)
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    pr(f"# l9_neargap [{variant}]: kinds {kinds}, Laplace order K{K} where computed; covered J = {sorted(mods)}; rows (Delta: K) = " +
       "; ".join(f"J={J}: " + ",".join(f"{D:.0f}:{k}" for D, k in sorted(mods[J][2].items())) for J in sorted(mods)))
    for Dg in Dgs:
        NG = NearGap(Dg); F = NG.F; pn, pw, Lm = NG.pn, NG.pw, NG.Lam
        blocks = {}
        for h, DF, r0, ref in load_states(Dg):
            mfl, J = r0['mfl'], r0['J']; u = 1/mfl**2; s = SIGMA[h]
            x = pn*np.sqrt(u); y = x*x
            P = remk.Pders_vec(Lm*pn, r0['eta'], nmax=0)[0]; jz = j0((2*J + 1)*x)
            pf2, wf2 = gl01(max(400, 2*J + 100))
            Js, B2, B3, _ = flat_parts_nf(F, h, u, J, pf2, wf2); i = int(np.searchsorted(Js, J)); B2, B3 = B2[i], B3[i]
            p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
            cv2 = (p2*2*(P - jz)).sum(); cv3 = (p3*u*s*(P - jz)).sum()
            fac = (1 - y)**2 if h == '++' else (1 - y)**-2
            rem = (p2*envk(2, x)*fac*r0['rho'][2]).sum() + (p3*u*s*envk(3, x)*fac*r0['rho'][3]).sum()
            Ns = abs(B2 + cv2) + abs(B3 + cv3)
            blocks.setdefault((round(mfl, 9), J, DF), {})[h] = dict(base=B2 + B3 + cv2 + cv3 + rem, Ns=Ns, r0=r0, x=x, y=y, p2=p2, p3=p3,
                                                                     u=u, m=r0['m'])
        recs = []
        for (mfl, J, DF), b in blocks.items():
            hs = [h for h in ('++', '+-') if h in b]
            base = np.array([b[h]['base'] for h in hs]); Ns = np.array([b[h]['Ns'] for h in hs])
            dF = np.zeros((len(hs), len(hs)))
            if J in mods and J % 2 == 0:
                e = b[hs[0]]; x, y, p2, p3, u, m = e['x'], e['y'], e['p2'], e['p3'], e['u'], e['m']
                fx = np.where(x < XC, TLOW[0], 1.0)
                G = mods[J][0](x, m)
                rho2, rho3 = e['r0']['rho'][2], e['r0']['rho'][3]
                r20 = np.polyfit(x[:3], rho2[:3], 2)[-1]
                for kk, h in enumerate(hs):
                    Dh = A.Dsh(h, J, x)
                    dF[kk, kk] = (fx*p2*A.env(h, x, 2)*(G[('C2', h)]/m**2 - (rho2 - Dh*r20))).sum()
                    if len(kinds) == 2:
                        dF[kk, kk] += (fx*p3*u*A.env(h, x, 3)*(G[(kinds[1], h)]/m**2 - A.SIG[h]*(rho3 - Dh*r20))).sum()
                if len(hs) == 2:
                    eo2 = np.sqrt(A.env('++', x, 2)*A.env('+-', x, 2))
                    o = (fx*p2*eo2*G[('C2', 'off')]).sum()/m
                    if len(kinds) == 2:
                        o += (fx*p3*u*eo2/(2 - y)*G[(kinds[1], 'off')]).sum()/m
                    dF[0, 1] = dF[1, 0] = o
            recs.append(dict(mfl=mfl, J=J, DF=DF, hs=hs, base=base, Ns=Ns, dF=dF))
        grid = sorted(set(r['mfl'] for r in recs))
        def lam(r, t):
            Fm = np.diag(r['base']) + t*r['dF']; n = 1/np.sqrt(r['Ns'])
            return np.linalg.eigvalsh(n[:, None]*Fm*n[None, :])[0]
        def nval(t):
            return next((n0 for n0 in grid if all(lam(r, t) > 0 for r in recs if r['mfl'] >= n0 - 1e-9)), None)
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
        pr(f"  Delta_g = {Dg} [{variant}]: n = {nval(1.0)} (no correction: {nval(0.0)}); min lam/N_s = {lam(w1, 1.0):+.4f} "
           f"[m/Lam={w1['mfl']} J={w1['J']} Df={w1['DF']} {w1['hs']}] (no corr.: {lam(w0, 0.0):+.4f} [m/Lam={w0['mfl']} J={w0['J']}]); "
           f"t* = ({tst[1]:+.1f}, {tst[0]:+.1f})")
        for J in sorted(mods):
            rj = [r for r in recs if r['J'] == J]
            if rj:
                wj = min(rj, key=lambda r: lam(r, 1.0))
                pr(f"      J={J:2d}: min lam/N_s {lam(wj, 1.0):+.4f} (uncorr. {min(lam(r, 0.0) for r in rj):+.4f}) at m/Lam={wj['mfl']}; "
                   f"max |dF_hh|/N_s {max(abs(r['dF'][k, k])/r['Ns'][k] for r in rj for k in range(len(r['hs']))):.3f}; "
                   f"max |dF_+-|/sqrt(NN) {max((abs(r['dF'][0, 1])/np.sqrt(r['Ns'].prod()) if len(r['hs']) == 2 else 0) for r in rj):.3f}")
        pr("      uncovered (unchanged) J: " + str(sorted(set(r['J'] for r in recs if r['J'] not in mods))))
    open(os.path.join(HERE, 'cache', 'l9_neargap_%s_%s.txt' % (variant, '_'.join(map(str, Dgs)))), 'w').write("\n".join(lines) + "\n")


if __name__ == '__main__':
    main(sys.argv[1], [int(a) for a in sys.argv[2:]] or [100, 300, 1000])
