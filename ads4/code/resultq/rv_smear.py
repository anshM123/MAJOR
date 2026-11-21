"""
rv_smear (RESULT_Q_REVISION.md §1): graviton-specific (beyond the scalar-skeleton + leading-Regge HYB model) corrections of the SMEARED
pieces of the Result Q functional Omega^Q (delta = 0.2, 'Q3-full-delta0.2'):  dB2 = int psi2 C'_EE (k = 2), dB3 = int psi3 C^_3 (k = 3).

Source of the correction: Result D's EE K2 crystal rows (the near-gap chain of EE_NEARGAP.md, run 'final_xc0' = convention A: q3/S1_S4
conventions, EE-own helicity basis, Hampel despike, c_off = m Gamma_off), read from spinfunc/cache/w2_recs_final_xc0.pkl (grids only,
read-only).  Per (h, J): Gamma_{C2,h}(x, m), Gamma_{C3H,h}(x, m) in c2 units x LO shape (l9_analysis.grid), cubic spline in x = nu/m,
linear in 1/m between Delta rows, held outside (l9_neargap.model, imported read-only).  The correction of a HYBRID state is (q3_analysis.
neargap with rho = 0):
    dF_hh = sum_p w psi2(p) env(h,x,2) Gamma_C2,h / m^2  +  sum_p w psi3(p) u env(h,x,3) Gamma_C3H,h / m^2,   x = p sqrt(u)
    dF_+- (J >= 4, off-diagonal) = sum w psi2 sqrt(env++ env+-) Gamma_C2,off/m + sum w psi3 u sqrt(env++ env+-)/(2-y) Gamma_C3H,off/m
with Q's psi2, psi3 in place of NF02's (same B2/B3 normalisation as stage6_model.hyb: B2 = int psi2 (2-y) D_h, B3 = sigma u int psi3 D_h).
Validation: the same code with NF02's psi's reproduces the stored Result D dF of the hybrid states (printed).

Variants (small x below the first row node, and the edge window x > x_last = 0.8):
   A  spline below / held above (Result D convention A)
   B  held below / held above
   Z  zero below x0 / zero above x_last (only the K2-covered window)
   E  as A, but above x_last the relative correction Gamma/D-shape is held (Gamma(x) = Gamma(x_last) env-shape ratio; +- only differs)
   X  as A, edge window x > 0.8 multiplied by 10 (stress)
Output: FF/N at the binding states with the measured forward s_d (BELOW_FOLD2/3, H0 gauge) and with dB2 + dB3 added; relative dB2/B2.
Usage: python3 rv_smear.py   -> cache/rv_smear.txt, cache/rv_smear.pkl
"""
import os, sys, io, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF); sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
import numpy as np
from scipy.interpolate import CubicSpline

CACHE = os.path.join(HERE, 'cache')
TAG = os.environ.get('RV_W2TAG', 'final_xc0')
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def env(h, x, k):
    y = x*x
    e = (2 - y) if h == '++' else (2 - y)/(1 - y)**4
    return e if k == 2 else e/(2 - y)


def Dh(h, J, x):
    import l7_dc2_analysis as An
    return An.Dh(h, J, x)


def load_models():
    import l9_neargap as NG
    g = pickle.load(open(os.path.join(SF, 'cache', 'w2_recs_%s.pkl' % TAG), 'rb'))['ee']['grids']
    em, om = {}, {}
    for J, rows in g['even'].items():
        rr = rows
        if J >= 4:
            rr = {D: {**r, **{(kd, 'off'): r[(kd, 'off')]*r['m'] for kd in ('C2', 'C3H')}} for D, r in rows.items()}
        if len(rr) == 1:
            D = list(rr)[0]; rr = {D: rr[D], D + 1e-6: dict(rr[D], m=rr[D]['m']*(1 + 1e-9))}
        em[J] = (NG.model(rr, J, ('C2', 'C3H'))[0], rr)
    for J, rows in g['odd'].items():
        Ds = sorted(rows); ms = np.array([rows[D]['m'] for D in Ds])
        spl = {D: {kd: CubicSpline(rows[D]['xs'], rows[D][kd], bc_type='natural') for kd in ('C2', 'C3H')} for D in Ds}
        xmax = {D: rows[D]['xs'][-1] for D in Ds}

        def f(x, m, Ds=Ds, ms=ms, spl=spl, xmax=xmax):
            if m <= ms[0]:
                i0 = i1 = 0; t = 0.0
            elif m >= ms[-1]:
                i0 = i1 = len(Ds) - 1; t = 0.0
            else:
                i1 = int(np.searchsorted(ms, m)); i0 = i1 - 1
                t = (1/ms[i0] - 1/m)/(1/ms[i0] - 1/ms[i1])
            return {kd: (1 - t)*spl[Ds[i0]][kd](np.minimum(x, xmax[Ds[i0]])) + t*spl[Ds[i1]][kd](np.minimum(x, xmax[Ds[i1]]))
                    for kd in ('C2', 'C3H')}
        om[J] = (f, rows)
    return em, om


def xrange_of(rows):
    return max(r['xs'][0] for r in rows.values()), min(r['xs'][-1] for r in rows.values())


def corr(model, rows, h, J, x, m, variant, odd=False):
    """Gamma's at the nodes x for the variant; returns dict kd -> array (diag of helicity h) and offs."""
    x0, x1 = xrange_of(rows)
    xe = np.clip(x, x0, None) if variant == 'B' else x
    G = model(xe, m)
    out = {}
    keys = ['C2', 'C3H'] if odd else [('C2', h), ('C3H', h)] + ([('C2', 'off'), ('C3H', 'off')] if J >= 4 and not odd else [])
    for k in keys:
        v = np.array(G[k], float)
        if not odd and k[1] == 'off':
            v = v/m                                                   # model holds c_off = m Gamma_off
        if variant == 'Z':
            v = np.where((x < x0) | (x > x1), 0.0, v)
        if variant == 'E' and (h == '+-'):
            xl = x1
            v = np.where(x > xl, v*(env(h, xl, 2)/env(h, np.minimum(x, 0.999), 2)), v)   # hold Gamma/env_+- (relative) beyond x_last
        if variant == 'X':
            v = np.where(x > x1, 10*v, v)
        out[k] = v
    return out


def dF_state(models, psi2, psi3, h, J, mfl, Lam, variant='A', n=256):
    em, om = models
    odd = J % 2 == 1
    if (odd and J not in om) or (not odd and J not in em):
        return None
    model, rows = (om[J] if odd else em[J])
    x_, w_ = np.polynomial.legendre.leggauss(n); pn = 0.5*(x_ + 1); pw = 0.5*w_
    u = 1/mfl**2; m = mfl*Lam
    x = pn*np.sqrt(u); y = x*x
    p2 = psi2(pn)*pw; p3 = psi3(pn)*pw
    G = corr(model, rows, h, J, x, m, variant, odd)
    if odd:
        d2 = (p2*env(h, x, 2)*G['C2']).sum()/m**2; d3 = (p3*u*env(h, x, 3)*G['C3H']).sum()/m**2
        return dict(d2=d2, d3=d3, off=0.0)
    d2 = (p2*env(h, x, 2)*G[('C2', h)]).sum()/m**2; d3 = (p3*u*env(h, x, 3)*G[('C3H', h)]).sum()/m**2
    off = 0.0
    if J >= 4:
        eo = np.sqrt(env('++', x, 2)*env('+-', x, 2))
        off = (p2*eo*G[('C2', 'off')]).sum()/m + (p3*u*eo/(2 - y)*G[('C3H', 'off')]).sum()/m
    return dict(d2=d2, d3=d3, off=off)


def validate(models):
    """NF02 psi's through the same code vs stored Result D dF of hybrid states (W2_GAUGE = none)."""
    from nf_common import NF
    F = NF()
    st = pickle.load(open(os.path.join(SF, 'cache', 'w2_recs_%s.pkl' % TAG), 'rb'))['ee']['recs']
    pr("## validation: NF02 psi's through rv_smear vs Result D's stored dF (hybrid states; Gauss 64 there, 256 here)")
    worst = 0.0
    for Dg in (300, 1000):
        Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
        for r in st[Dg]:
            if r['DF'] != 'hyb' or r['J'] not in (0, 2, 4, 5) or r['mfl'] > 3.01:
                continue
            for kk, h in enumerate(r['hs']):
                d = dF_state(models, F.psi2, F.psi3, h, r['J'], r['mfl'], Lam, 'A', n=64)
                mine = d['d2'] + d['d3']
                ref = r['dF'][kk, kk]
                worst = max(worst, abs(mine - ref)/max(abs(ref), 1e-12))
                pr(f"   Dg={Dg} J={r['J']} m/Lam={r['mfl']} {h}: rv {mine:+.6e}  stored {ref:+.6e}" +
                   (f"  | off rv {d['off']:+.3e} stored {r['dF'][0, 1]:+.3e}" if len(r['hs']) == 2 and kk == 0 else ''))
    pr(f"   max relative deviation (diag): {worst:.1e}")


def main():
    import stage6_model as M, common_q as CQ
    from stage6_scan import BMAX
    models = load_models()
    validate(models)
    F, c2, c3 = M.load('Q3-full-delta0.2'); a, c = F['a'], F['c']
    psi2 = lambda p: np.polynomial.polynomial.polyval(p, c2)
    psi3 = lambda p: np.polynomial.polynomial.polyval(p, c3)
    # measured forward s_d (H0 gauge; BELOW_FOLD2/3): ++ J = 0, 2, 4: +0.509/+0.509/+0.517 (300), +0.508/+0.508/+0.510 (1000); +- J = 4: -0.516/-0.510
    SD = {(300, '++', 0): 0.509, (300, '++', 2): 0.509, (300, '++', 4): 0.517, (300, '+-', 4): -0.516,
          (1000, '++', 0): 0.508, (1000, '++', 2): 0.508, (1000, '++', 4): 0.510, (1000, '+-', 4): -0.510}
    out = {}
    pr("\n## Result Q (delta = 0.2) binding states: HYB parts, measured forward s_d, and the K2-chain graviton correction of the smeared parts")
    pr("   columns: B2, B3 (HYB) | N | FF/N: HYB K=0 -> + measured s_d -> + dB2 + dB3 (variant A) | dB2/B2, dB3/|B3| | variants A/B/Z/E/X of (dB2+dB3)/N")
    for Dg in (300, 1000):
        Lam = M.lam(Dg); tab = M.Tables('Q3-full-delta0.2', Lam, BMAX[Dg])
        mg = lambda J: np.sqrt(1/M.u_gap(J, Dg))
        states = [('++', 0, 1.0 + 1e-12), ('++', 0, 1.05), ('++', 2, mg(2)), ('++', 2, 1.05), ('++', 2, 1.118), ('++', 4, mg(4)),
                  ('+-', 4, mg(4)), ('+-', 4, 1.05), ('+-', 5, mg(5)), ('+-', 5, 1.15), ('+-', 5, 1.2)]
        for h, J, mfl in states:
            u = 1/mfl**2
            Js, B2, B3 = CQ.parts(h, u, max(J, 4), [c2], [c3]); k = list(Js).index(J)
            FF, N, _, parts = M.hyb(F, h, u, float(J), B2[k, 0], B3[k, 0], tab, Lam)
            sd = SD.get((Dg, h, J), 0.0)
            FFs = FF - a*2*u**3*sd
            res = {v: dF_state(models, psi2, psi3, h, J, mfl, Lam, v) for v in ('A', 'B', 'Z', 'E', 'X')}
            if res['A'] is None:
                pr(f"   Dg={Dg} {h} J={J} m/Lam={mfl:.4f}: no K2 grid"); continue
            dA = res['A']['d2'] + res['A']['d3']
            B2h, B3h = float(parts['B2']), float(parts['B3'])
            Nn = abs(B2h + res['A']['d2']) + abs(B3h + res['A']['d3']) + abs(float(parts['FW'])) + c*abs(float(parts['B4']))
            pr(f"   Dg={Dg} {h} J={J} m/Lam={mfl:.4f}: B2 {B2h:+.4f} B3 {B3h:+.4f} | N {float(N):.4f} | FF/N {float(FF/N):+.4f} -> "
               f"{float(FFs/N):+.4f} (s_d {sd:+.3f}) -> {float((FFs + dA)/Nn):+.4f} | dB2/B2 {res['A']['d2']/B2h:+.2e} "
               f"dB3/|B3| {res['A']['d3']/abs(B3h):+.2e} | (dB2+dB3)/N: " +
               " ".join(f"{v}:{(res[v]['d2'] + res[v]['d3'])/float(N):+.2e}" for v in res) +
               (f" | off/N {res['A']['off']/float(N):+.1e}" if J >= 4 and J % 2 == 0 else ''))
            out[(Dg, h, J, round(mfl, 4))] = dict(B2=B2h, B3=B3h, N=float(N), FF=float(FF), FFs=float(FFs), res=res)
    open(os.path.join(CACHE, 'rv_smear.txt'), 'w').write(buf.getvalue())
    pickle.dump(out, open(os.path.join(CACHE, 'rv_smear.pkl'), 'wb'))


if __name__ == '__main__':
    main()
