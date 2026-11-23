"""
p2spin/s4_jstar.py (P2_SPINNING.md §5): Theorem W's J* with CLASS-WISE premise constants.
Same criterion and the same read-only cost table as lsw/lsw_jstar.py (lsw/cache/costs_1000.json):
    sum_s Lam_s |kappa_s^sc| c_s(form)  +  K_sh sum c_sh(R2)  +  K4 sum_{4a,4b} c_4(R2)/m^2  <=  theta * H.
Lam_s per class s in {A, B, C1, C2, C3} (multiples of CHMRS's scalar kappa_s^sc).  The lsw_jstar rows are reproduced first (check).
Usage: python3 s4_jstar.py  -> cache/s4_jstar.txt
"""
import json, os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, '..', 'lsw', 'cache', 'costs_1000.json')))
SC = {'A ': 86.75/6 - 8.0, 'B ': 11/6, 'C1': 1/6, 'C2': 1/6, 'C3': 11.75/6}
out = []; pr = lambda s: (print(s, flush=True), out.append(s))


def load(r, form, lam):
    c = r['cost']
    return (sum(lam[k]*abs(SC[k])*abs(c[k][form]) for k in SC),
            sum(abs(c[k][1]) for k in c if k.startswith('sh')), sum(abs(c[k][1]) for k in ('4a', '4b'))/r['m']**2)


def H(r):
    f = {h: v['floor'] for h, v in r['fl'].items()}
    return f['++']*f['+-']/(f['++'] + f['+-']) if len(f) == 2 else list(f.values())[0]


Js = sorted(set(r['J'] for r in R))


def jstar(form, lam, Ksh, K4, th):
    fails = sorted(set(r['J'] for r in R if sum(w*x for w, x in zip((1, Ksh, K4), load(r, form, lam))) > th*H(r)))
    ok = [J for J in Js if all(J2 < J or J2 not in fails for J2 in Js)]
    return (min(ok) if ok else 'none'), fails


def uni(L):
    return {k: L for k in SC}


if __name__ == '__main__':
    pr('# check against lsw/cache/jstar_1000.txt (uniform Lam_amp):')
    for form, lab in ((0, 'R1'), (1, 'R2')):
        for La, Ksh, K4 in ((1, 10, 1), (3, 10, 10)):
            for th in (1.0, 0.5):
                js, f = jstar(form, uni(La), Ksh, K4, th)
                pr(f'  {lab} Lam={La} K_sh={Ksh} K4={K4} theta={th}: J* = {js}  fails {f}')
    # model premise constants of P2_SPINNING.md §3.2 (direct / conj pairing of the helicity-diagonal kernel model [M])
    cfg = [{"label": "model-direct, K_sh=10, K4=1", "lam": {"A": 3.17, "B": 1.0, "C1": 1.0, "C2": 1.0, "C3": 2.03}, "Ksh": 10, "K4": 1}, {"label": "model-direct, K_sh=10, K4=10", "lam": {"A": 3.17, "B": 1.0, "C1": 1.0, "C2": 1.0, "C3": 2.03}, "Ksh": 10, "K4": 10}, {"label": "model-direct, K_sh=30, K4=10", "lam": {"A": 3.17, "B": 1.0, "C1": 1.0, "C2": 1.0, "C3": 2.03}, "Ksh": 30, "K4": 10}, {"label": "model-conj, K_sh=10, K4=10", "lam": {"A": 11.5, "B": 1.0, "C1": 205.0, "C2": 1.0, "C3": 6.4}, "Ksh": 10, "K4": 10}, {"label": "model-direct, K_sh=50, K4=10", "lam": {"A": 3.17, "B": 1.0, "C1": 1.0, "C2": 1.0, "C3": 2.03}, "Ksh": 50, "K4": 10}]
    for c in cfg:
        lam = {k: c['lam'].get(k.strip(), c['lam'].get(k, 1.0)) for k in SC}
        for form, lab in ((0, 'R1 (de-shift)'), (1, 'R2 (stripped)')):
            for th in (1.0, 0.5):
                js, f = jstar(form, lam, c['Ksh'], c['K4'], th)
                pr(f'  [{c["label"]}] {lab} Lam={ {k.strip(): round(v, 3) for k, v in lam.items()} } K_sh={c["Ksh"]} K4={c["K4"]} '
                   f'theta={th}: J* = {js}  fails {f}')
    open(os.path.join(HERE, 'cache', 's4_jstar.txt'), 'w').write('\n'.join(out) + '\n')
