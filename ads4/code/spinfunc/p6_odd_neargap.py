"""
p6 (PARITY_ODD_ROUTEA.md sect. 6): route-A' near-gap positivity of the ODD-J (+-, parity-odd exchange) states with the computed
graviton-specific corrections of Ĉ_2 and Ĉ_3 (p4 grids), in exactly the conventions of l9_neargap (variant C23H_K2):
  F = base_{+-}  (stripped flat B2 + B3, (P - J0) terms, exact scalar remainder x envelopes; no inversion; = h2/l8/l9 base)
    + sum_nodes psi2 w env_{2,+-} [Gamma_{2}(x; m)/m^2 - (rho2 - D_+- rho2(0))]
    + sum_nodes psi3 w u env_{3,+-} [Gamma_{3}(x; m)/m^2 + (rho3 - D_+- rho2(0))]       (sigma_{+-} = -1)
  (1x1: no helicity mixing for odd J -- the odd structure feeds +- only, and EO interference cancels in the S^0 sum, p3).
Direct values on the computed Delta rows, linear in 1/m between rows, held beyond; x beyond the last computed x held.
t*: scale factor on the whole correction at which the odd sector first fails n = 1 (lambda < 0 for some state with m/Lam >= grid min).
Usage: python3 p6_odd_neargap.py Dg [Dg ...]   (P4_DESHIFT=0: literal CLS-D, no de-shift)
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
import p4_odd_analysis as A

DMAX = float(os.environ.get('P6_DMAX', '3200'))
KINDS = tuple(os.environ.get('P6_KINDS', 'C2,C3H').split(','))
K3ROW = os.environ.get('P6_K3', '0') == '1'      # use the K3 grid for the Delta = 100 row where computed


def model(J):
    g = pickle.load(open(os.path.join(HERE, 'cache', 'p4_grid%s.pkl' % ('' if A.DESHIFT else '_noshift')), 'rb'))[J]
    if K3ROW:
        g3 = pickle.load(open(os.path.join(HERE, 'cache', 'p4_grid_k3%s.pkl' % ('' if A.DESHIFT else '_noshift')), 'rb')).get(J, {})
        g = dict(g); g.update(g3)
    Ds = sorted(D for D in g if D <= DMAX); ms = np.array([g[D]['m'] for D in Ds])
    spl = {D: {kd: CubicSpline(g[D]['xs'], g[D][kd], bc_type='natural') for kd in KINDS} for D in Ds}
    xmax = {D: g[D]['xs'][-1] for D in Ds}
    def f(x, m):
        if m <= ms[0]:
            i0 = i1 = 0; t = 0.0
        elif m >= ms[-1]:
            i0 = i1 = len(Ds) - 1; t = 0.0
        else:
            i1 = int(np.searchsorted(ms, m)); i0 = i1 - 1
            t = (1/ms[i0] - 1/m)/(1/ms[i0] - 1/ms[i1])
        return {kd: (1 - t)*spl[Ds[i0]][kd](np.minimum(x, xmax[Ds[i0]])) + t*spl[Ds[i1]][kd](np.minimum(x, xmax[Ds[i1]]))
                for kd in KINDS}
    return f, Ds


def main(Dgs, Js=(5, 7, 9)):
    mods = {}
    for J in Js:
        try:
            mods[J] = model(J)
        except Exception as ex:
            print(f"J={J}: no grid ({ex})")
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    pr(f"# p6 [kinds {KINDS}, K3 row at Delta=100: {K3ROW}]: odd-J (+-, parity-odd) states, corrections, K2 linear chain, de-shift {A.DESHIFT}; covered J = "
       + ", ".join(f"{J} (rows {','.join(f'{D:.0f}' for D in mods[J][1])})" for J in sorted(mods)))
    for Dg in Dgs:
        NG = NearGap(Dg); F = NG.F; pn, pw, Lm = NG.pn, NG.pw, NG.Lam
        recs = []
        allmin = []
        for h, DF, r0, ref in load_states(Dg):
            mfl, J = r0['mfl'], r0['J']
            u = 1/mfl**2; s = SIGMA[h]
            x = pn*np.sqrt(u); y = x*x
            P = remk.Pders_vec(Lm*pn, r0['eta'], nmax=0)[0]; jz = j0((2*J + 1)*x)
            pf2, wf2 = gl01(max(400, 2*J + 100))
            Js_, B2, B3, _ = flat_parts_nf(F, h, u, J, pf2, wf2); i = int(np.searchsorted(Js_, J)); B2, B3 = B2[i], B3[i]
            p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
            cv2 = (p2*2*(P - jz)).sum(); cv3 = (p3*u*s*(P - jz)).sum()
            fac = (1 - y)**2 if h == '++' else (1 - y)**-2
            rem = (p2*envk(2, x)*fac*r0['rho'][2]).sum() + (p3*u*s*envk(3, x)*fac*r0['rho'][3]).sum()
            Ns = abs(B2 + cv2) + abs(B3 + cv3)
            base = B2 + B3 + cv2 + cv3 + rem
            allmin.append((base/Ns, mfl, J, h, DF))
            if J % 2 == 0:
                continue
            dF = 0.0
            if J in mods:
                m = r0['m']
                G = mods[J][0](x, m)
                rho2, rho3 = r0['rho'][2], r0['rho'][3]
                r20 = np.polyfit(x[:3], rho2[:3], 2)[-1]
                Dh = A.Dc('+-', J, x)
                dF = (p2*A.env('+-', x, 2)*(G['C2']/m**2 - (rho2 - Dh*r20))).sum()
                if 'C3H' in KINDS:
                    dF += (p3*u*A.env('+-', x, 3)*(G['C3H']/m**2 - A.SIG['+-']*(rho3 - Dh*r20))).sum()
            recs.append(dict(mfl=mfl, J=J, DF=DF, base=base, Ns=Ns, dF=dF, cov=J in mods))
        grid = sorted(set(r['mfl'] for r in recs))
        lam = lambda r, t: (r['base'] + t*r['dF'])/r['Ns']
        def ok(t):
            return all(lam(r, t) > 0 for r in recs if r['mfl'] >= grid[0] - 1e-9)
        w0 = min(recs, key=lambda r: lam(r, 0.0)); w1 = min(recs, key=lambda r: lam(r, 1.0))
        tst = []
        for sgn in (1, -1):
            lo, hi = 0.0, 1.0
            while ok(sgn*hi) and hi < 1e4:
                lo, hi = hi, 2*hi
            for _ in range(40):
                mid = 0.5*(lo + hi)
                lo, hi = (mid, hi) if ok(sgn*mid) else (lo, mid)
            tst.append(sgn*lo)
        gm = min(allmin)
        pr(f"  Delta_g = {Dg}: odd-J states {len(recs)} (covered {sum(r['cov'] for r in recs)}); odd sector positive with correction: {ok(1.0)}; "
           f"min F/N_s = {lam(w1, 1.0):+.4f} [m/Lam={w1['mfl']} J={w1['J']} Df={w1['DF']}] (uncorrected {lam(w0, 0.0):+.4f} "
           f"[m/Lam={w0['mfl']} J={w0['J']}]); t* = ({tst[1]:+.1f}, {tst[0]:+.1f}); all-state uncorrected min {gm[0]:+.4f} {gm[1:]}")
        for J in sorted(set(r['J'] for r in recs)):
            rj = [r for r in recs if r['J'] == J]
            wj = min(rj, key=lambda r: lam(r, 1.0))
            pr(f"      J={J:2d}{'' if J in mods else ' (uncovered)'}: min F/N_s {lam(wj, 1.0):+.4f} (uncorr. {min(lam(r, 0.0) for r in rj):+.4f}) "
               f"at m/Lam={wj['mfl']} Df={wj['DF']}; max |dF|/N_s {max(abs(r['dF'])/r['Ns'] for r in rj):.4f}; "
               f"binding-state dF/N_s {wj['dF']/wj['Ns']:+.4f}")
    open(os.path.join(HERE, 'cache', 'p6_odd_neargap%s_%s%s_%s.txt' % ('' if A.DESHIFT else '_noshift', '_'.join(KINDS), '_k3' if K3ROW else '',
                                                                    '_'.join(map(str, Dgs)))), 'w').write(
        "\n".join(lines) + "\n")


if __name__ == '__main__':
    main([int(a) for a in sys.argv[1:]] or [300, 1000, 100])
