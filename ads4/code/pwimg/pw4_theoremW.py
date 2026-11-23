"""
pwimg/pw4_theoremW.py (PW_IMAGES.md §5): Theorem W without premise P3 (R1).  By pw1 a stripped (R2) raw correction cannot occur:
every raw correction is A_g[F] with F even and strip-analytic.  The only way to leave the EE de-shift form is an F whose ratio to c
does not vanish at the Regge-spin-0 point mu = +-i/2.  The O(1/m^2) mechanisms (radial-weight shifts of the leading kernel by
|x|^sA |y|^sB, pw2) produce F = c_s(mu) S(mu; beta) with the EE/scalar ratio c_s of pw2 instead of c.  Exact cost (pw3):
      R1'(S, s) = 2 [ sum psi2 (c_s/c) S dnu/Lam + D c_s(i/2) S(i/2) ],    S(i/2): A -> eta, B -> -eta/4, C1, C2, C3, 4b -> 0,
                  4a -> eta^2/16, sh+-1 -> (1 + eta)/2, sh+-2 -> (eta + P_2(eta))/2   (P_{-1/2+i mu}(cosh b) = 1 at mu = +-i/2).
'R1'' form' := max over s in SHIFTS of |R1'(S, s)|  (sign-free; s = (0,0) is the EE de-shift form R1 itself).
Then Theorem W's criterion (lsw_theorem/lsw_jstar, read-only input cache lsw/cache/costs_1000.json for states and floors) is re-run
with R1' in place of R2 for every class.  Output: cache/pw4_theoremW.txt, cache/pw4_theoremW.json
"""
import os, sys, io, json
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'lsw'))
import numpy as np
import lsw_common as W          # read-only (its own read-only imports of n2sub, largespin)
from pw_common import cfac

OUT = io.StringIO()
pr = lambda s: (print(s, flush=True), OUT.write(s + '\n'))
D = json.load(open(os.path.join(HERE, 'cache', 'pw3_cost.json')))['1000']['D']
# exact EE/scalar ratios for radial shifts (sA, sB) (pw2_radial.py, cache/pw2_radial.txt)
CS = {
    (0, 0): lambda v: cfac(v),
    (1, 1): lambda v: 32*(16*v**4 + 120*v**2 + 119)/((4*v**2 + 49)*(4*v**2 + 121)),
    (0, 2): lambda v: 32*(4*v**2 + 9)*(4*v**2 + 49)/((4*v**2 + 81)*(4*v**2 + 169)),
    (2, 2): lambda v: 32*(4*v**2 + 9)*(4*v**2 + 49)/((4*v**2 + 81)*(4*v**2 + 169)),
    (3, 3): lambda v: 32*(16*v**4 + 376*v**2 + 1143)/((4*v**2 + 121)*(4*v**2 + 225)),
}
SHIFTS = [(1, 1), (0, 2), (3, 3)]
SC = {'A ': 86.75/6 - 8.0, 'B ': 11/6, 'C1': 1/6, 'C2': 1/6, 'C3': 11.75/6}
R = json.load(open(os.path.join(W.CACHE, 'costs_1000.json')))


def S_at_half(name, beta):
    eta = np.cosh(beta)
    return {'A ': eta, 'B ': -eta/4, 'C1': 0.0, 'C2': 0.0, 'C3': 0.0, '4a': eta**2/16, '4b': 0.0,
            'sh+1': (1 + eta)/2, 'sh-1': (1 + eta)/2, 'sh+2': (eta + (3*eta**2 - 1)/2)/2, 'sh-2': (eta + (3*eta**2 - 1)/2)/2}[name]


def r1p(Sv, Sh, s):
    ratio = CS[s](W.NU)/cfac(W.NU)
    return 2*((W.PS*ratio*Sv).sum()*W.DN/W.L + D*CS[s](0.5j).real*Sh)


def main():
    pr(__doc__.split('Output')[0])
    pr('c_s(i/2) (monopole weight per unit S(i/2)): ' + ', '.join(f'{s}: {CS[s](0.5j).real:+.4f}' for s in CS))
    rows = []
    for r in R:
        b = r['beta']; m = r['m']
        S = W.structures(b)
        c = {}
        for k, Sv in S.items():
            Sh = S_at_half(k, b)
            vals = [r1p(Sv, Sh, s) for s in SHIFTS]
            mono = [2*D*CS[s](0.5j).real*Sh for s in SHIFTS]
            c[k] = dict(R1=abs(r['cost'][k][0]), R2=abs(r['cost'][k][1]), R1p=max(abs(v) for v in vals)/m**2,
                        mono=max(abs(v) for v in mono)/m**2)
        f = {h: v['floor'] for h, v in r['fl'].items()}
        H = f['++']*f['+-']/(f['++'] + f['+-']) if len(f) == 2 else list(f.values())[0]
        rows.append(dict(J=r['J'], mfl=r['mfl'], m=m, beta=b, H=H, c=c))
    pr('\nper-unit costs relative to H (= fraction of the 2x2 budget per unit c2-coefficient), min-H mass per J shown:')
    pr('  J    m/L    beta  | class: R1 (EE de-shift) / R2 (stripped, now excluded) / R1\' (c-ratio, exact, incl. monopole) / monopole part')
    byJ = {}
    for q in rows:
        byJ.setdefault(q['J'], []).append(q)
    for J in sorted(byJ, key=lambda j: (j % 2, j)):
        q = min(byJ[J], key=lambda z: z['H'])
        pr(f"  {J:3d}  {q['mfl']:.3f}  {q['beta']:.3f} | " + '  '.join(
            f"{k.strip()}: {q['c'][k]['R1']/q['H']:.1e}/{q['c'][k]['R2']/q['H']:.1e}/{q['c'][k]['R1p']/q['H']:.1e}/{q['c'][k]['mono']/q['H']:.0e}"
            for k in ('A ', 'B ', 'C2', 'sh+1', 'sh+2')))
    # class-wise tolerances (min over masses), as LARGESPIN_W §4.2 but with R1'
    pr('\nclass-wise K_tol = H / cost (min over sampled masses):   J | A | B | C1 | C2 | C3 | index-shifted (sum) | scalar vector x1')
    tol = {}
    for J in sorted(byJ, key=lambda j: (j % 2, j)):
        t = {}
        for form in ('R1', 'R2', 'R1p'):
            t[form] = {k: min(q['H']/max(q['c'][k][form], 1e-300) for q in byJ[J]) for k in SC}
            t[form]['sh'] = min(q['H']/sum(q['c'][k][form] for k in q['c'] if k.startswith('sh')) for q in byJ[J])
            t[form]['vec'] = min(q['H']/sum(abs(SC[k])*q['c'][k][form] for k in SC) for q in byJ[J])
        tol[J] = t
        pr(f"  {J:3d} | " + ' | '.join(f"{k.strip()}: {t['R1'][k]:.3g}/{t['R2'][k]:.3g}/{t['R1p'][k]:.3g}" for k in ('A ', 'B ', 'C1', 'C2', 'C3', 'sh', 'vec')))
    pr('  (entries R1 / R2 / R1\')')
    # J* (lsw_jstar logic), with form R1 (P3 assumed), R2 (P3 dropped, old worst case), R1' (P3 dropped, exact worst case)
    pr('\nJ* (smallest sampled J0 with every sampled state J >= J0 passing; Lam_amp sum|kappa^sc| c + K_sh sum c_sh + K4 sum c_4/m^2 <= theta H):')
    Js = sorted(byJ)
    jres = {}
    for form, lab in (('R1', 'P3 assumed (EE de-shift, as W §4.3)'), ('R2', 'P3 dropped, stripped worst case (W §9.3)'),
                      ('R1p', 'P3 dropped, exact c-ratio worst case (this note)')):
        for (La, Ksh, K4) in ((1, 10, 1), (3, 10, 10), (10, 30, 10), (30, 100, 30)):
            for th in (1.0, 0.5):
                fails = []
                for q in rows:
                    shform = 'R2' if form == 'R2' else form
                    load = La*sum(abs(SC[k])*q['c'][k][form] for k in SC) + Ksh*sum(q['c'][k][shform] for k in q['c'] if k.startswith('sh')) \
                        + K4*sum(q['c'][k]['R2' if form == 'R2' else form] for k in ('4a', '4b'))/q['m']**2
                    if load > th*q['H']:
                        fails.append(q['J'])
                fails = sorted(set(fails))
                ok = [J for J in Js if all(J2 < J or J2 not in fails for J2 in Js)]
                js = min(ok) if ok else None
                jres[f'{form},{La},{Ksh},{K4},{th}'] = dict(Jstar=js, fails=fails)
                pr(f'  {lab:52s} (Lam_amp, K_sh, K4) = ({La}, {Ksh}, {K4}), theta = {th}: J* = {js}; failing sampled J: {fails}')
    json.dump(dict(rows=rows, tol={str(k): v for k, v in tol.items()}, jstar=jres), open(os.path.join(HERE, 'cache', 'pw4_theoremW.json'), 'w'))
    open(os.path.join(HERE, 'cache', 'pw4_theoremW.txt'), 'w').write(OUT.getvalue())


if __name__ == '__main__':
    main()
