"""
p2con/pc_jstar.py (P2_CONTRACT.md §4): Theorem W's J* with the COMPUTED O(1/m^2) block-side corrections (pc_run) in place of the
model constants Lam_amp * |kappa^sc|.

Per sampled state (same states, floors and per-unit costs as LARGESPIN_W / PW_IMAGES; read-only: lsw/cache/costs_1000.json,
pwimg/cache/pw4_theoremW.json):
  d2  = max over the entries (++,++), (+-,+-), (++,+-) [even J] or (+-,+-) [odd J] of |exact cost of the O(1/m^2) block-side correction|
        (exact R1'-type cost incl. the monopole; pc_run)                                               [computed, DV]
  dsh = K_sh * sum_{sh +-1, +-2} c_sh(R1')                                                             [premise P2-sh, M]
  d4  = K4 * sum_{4a, 4b} c_4(R1') / m^2                                                               [premise P4]
Criterion (Theorem W, sign-free): d2 + dsh + d4 <= theta H,  H = f++ f+- / (f++ + f+-) (even J) or f+- (odd J).
Variant 'signed 2x2': diag(f) + (C_exact - (dsh + d4))/theta PSD with |offdiag| + (dsh + d4) (uses the computed signs).
Also: Lam_eff = d2 / d2^sc with d2^sc = sum_s |kappa_s^sc| c_s(R1') (Theorem W's Lam_amp = 1 load), per state.
Map: 'same' = +-, 'opp' = ++ (P2_SPINNING §2.2).
Output: cache/jstar.txt, cache/jstar.json
"""
import os, sys, json, glob, io
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
from scipy.special import j0, j1
sys.path.insert(0, HERE)
import pc_run as Rn

OUT = io.StringIO()
pr = lambda s: (print(s, flush=True), OUT.write(s + '\n'))
SC = {'A ': 86.75/6 - 8.0, 'B ': 11/6, 'C1': 1/6, 'C2': 1/6, 'C3': 11.75/6}


def load_runs(pat='run_*_*.json'):
    rows = []
    for fn in sorted(glob.glob(os.path.join(HERE, 'cache', pat))):
        rows += json.load(open(fn))
    return {(r['J'], round(r['mfl'], 6)): r for r in rows}


C1PM = 31.8339     # flat O(y) J1/b coefficient of the +- partial wave (pc_flat / pc_p4 fit; = -kC(+-) of pc_run at beta -> 0)


def flat_c1(h, bb):
    """O(y) flat-tower function of the HYBRID floor's flat part (stripped, as coded in omegau flat_parts_nf), per helicity."""
    if h == '++':
        return -j0(bb)/6 - bb*j1(bb)/6 - j1(bb)/(6*bb)
    return (23/6)*j0(bb) - bb*j1(bb)/6 + C1PM*j1(bb)/bb


def bh_cost(xs, beta, m, J, mfl, h):
    """beyond-hybrid cost: 2[sum psi2 (f_total - nu^2 c1_flat(bhat p)) + D Q_total]/m^2 (f_flat entire: no defect)."""
    nu = Rn.W.NU
    f = 0; Q = 0
    for x in xs:
        F = Rn.Feval(np.array(x), beta, nu)[0]
        f = f + F/Rn.cfac(nu); Q += Rn.F_half(np.array(x), beta)[0].real
    if h is not None:
        bb = (2*J + 1)/mfl*nu/Rn.W.L
        f = f - nu**2*flat_c1(h, bb)
    return 2*((Rn.W.PS*f).real.sum()*Rn.W.DN/Rn.W.L + Rn.D1000*Q)/m**2


def main():
    pw4 = json.load(open(os.path.join(HERE, '..', 'pwimg', 'cache', 'pw4_theoremW.json')))['rows']
    runs = load_runs(); ksh = load_runs('ksh_*_*.json')
    states = []
    for q in pw4:
        key = (q['J'], round(q['mfl'], 6))
        if key not in runs:
            continue
        r = runs[key]; m = q['m']; c = q['c']
        e = r['ent']; kk = ksh[key]['ent']
        if q['J'] % 2 == 0:
            pairs = {('+-', '+-'): ('same|same', '+-'), ('++', '++'): ('opp|opp', '++'), ('++', '+-'): ('same|opp', None)}
        else:
            pairs = {('+-', '+-'): ('odd|odd', '+-')}
        Cblk = {k: e[v + '|full']['total'] for k, (v, h) in pairs.items()}
        Cker = {k: kk[v]['total'] for k, (v, h) in pairs.items()}
        C = {k: Cblk[k] + Cker[k] for k in pairs}
        Cbh = {k: bh_cost([e[v + '|full']['x'], kk[v]['x']], r['beta'], m, q['J'], q['mfl'], h) for k, (v, h) in pairs.items()}
        Ctot_chk = {k: bh_cost([e[v + '|full']['x'], kk[v]['x']], r['beta'], m, q['J'], q['mfl'], None) for k, (v, h) in pairs.items()}
        mode = os.environ.get('PC_MODE', 'full')
        Cuse = {'full': C, 'block': Cblk, 'bh': Cbh}[mode]
        d2 = max(abs(v) for v in Cuse.values())
        d2sc = sum(abs(SC[k])*c[k]['R1p'] for k in SC)
        csh = sum(c[k]['R1p'] for k in c if k.startswith('sh'))
        c4 = sum(c[k]['R1p'] for k in ('4a', '4b'))/m**2
        f = r['fl']
        states.append(dict(J=q['J'], mfl=q['mfl'], m=m, beta=q['beta'], H=q['H'], f=f, C={'%s,%s' % k: v for k, v in Cuse.items()}, Cfull={'%s,%s' % k: v for k, v in C.items()},
                           Cblk={'%s,%s' % k: v for k, v in Cblk.items()}, Cker={'%s,%s' % k: v for k, v in Cker.items()},
                           Cbh={'%s,%s' % k: v for k, v in Cbh.items()}, chk=max(abs(Ctot_chk[k] - C[k])/max(abs(C[k]), 1e-300) for k in C), d2=d2, d2sc=d2sc,
                           csh=csh, c4=c4, Leff=d2/d2sc))
    pr(__doc__.split('Output')[0])
    pr(f"{len(states)} states; mode = {os.environ.get('PC_MODE', 'full')} (full: block + kernel O(r), Theorem W convention; block: block side only; "
       f"bh: beyond-hybrid, flat O(y) tower of the hybrid floor subtracted); max consistency |recomputed/total - 1| = {max(s['chk'] for s in states):.1e}\n")
    pr('  J   m/L   beta  |  H        | d2 (computed)/H | d2^sc(Lam=1)/H | Lam_eff = d2/d2^sc | entries C/H: (+-,+-) (++,++) (++,+-) | K_sh tol | K4 tol (alone)')
    for s in sorted(states, key=lambda z: (z['J'] % 2, z['J'], z['mfl'])):
        Cs = ' '.join(f"{v/s['H']:+.2e}" for v in s['C'].values())
        s['Ksh_tol'] = (s['H'] - s['d2'])/s['csh']; s['K4_tol'] = (s['H'] - s['d2'])/s['c4']
        pr(f"  {s['J']:3d} {s['mfl']:.3f} {s['beta']:.4f} | {s['H']:.3e} | {s['d2']/s['H']:.3e} | {s['d2sc']/s['H']:.3e} | {s['Leff']:8.3f} | {Cs} | "
           f"{s['Ksh_tol']:9.3g} | {s['K4_tol']:9.3g}")
    Js = sorted(set(s['J'] for s in states))

    def jstar(Ksh, K4, th, signed=False):
        fails = []
        for s in states:
            dp = Ksh*s['csh'] + K4*s['c4']
            if not signed:
                ok = s['d2'] + dp <= th*s['H']
            else:
                f = s['f']
                if len(f) == 2:
                    a = f['++'] + (s['C']['++,++'] - dp)/th; cc = f['+-'] + (s['C']['+-,+-'] - dp)/th; b = (abs(s['C']['++,+-']) + dp)/th
                    ok = a > 0 and cc > 0 and a*cc >= b*b
                else:
                    ok = list(f.values())[0] + (s['C']['+-,+-'] - dp)/th > 0
            if not ok:
                fails.append(s['J'])
        fails = sorted(set(fails))
        good = [J for J in Js if all(J2 < J or J2 not in fails for J2 in Js)]
        return (min(good) if good else None), fails
    pr('\nJ* (smallest sampled J0 such that every sampled state with J >= J0 passes):')
    res = {}
    for signed in (False, True):
        for (Ksh, K4) in ((0, 0), (0, 1), (0, 4), (0, 10), (0, 35), (10, 4), (10, 10), (30, 10)):
            for th in (1.0, 0.5):
                js, fl = jstar(Ksh, K4, th, signed)
                res[f'{signed},{Ksh},{K4},{th}'] = dict(Jstar=js, fails=fl)
                pr(f"  {'signed 2x2' if signed else 'sign-free  '} K_sh = {Ksh:4g}, K4 = {K4:4g}, theta = {th}: J* = {js}; failing sampled J: {fl}")
    pr('\nLam_eff over states: min %.3f, max %.3f (even J), %s (odd J)' % (
        min(s['Leff'] for s in states if s['J'] % 2 == 0), max(s['Leff'] for s in states if s['J'] % 2 == 0),
        ', '.join('%.3f' % s['Leff'] for s in states if s['J'] % 2)))
    json.dump(dict(states=states, jstar=res), open(os.path.join(HERE, 'cache', 'jstar_%s.json' % os.environ.get('PC_MODE', 'full')), 'w'))
    open(os.path.join(HERE, 'cache', 'jstar_%s.txt' % os.environ.get('PC_MODE', 'full')), 'w').write(OUT.getvalue())


if __name__ == '__main__':
    main()
