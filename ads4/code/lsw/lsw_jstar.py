"""
lsw/lsw_jstar.py (LARGESPIN_W.md §4.3): J* of Theorem W for stated premise constants.
State passes iff  Lam_amp * sum_s |kappa_s^sc| c_s(form)  +  K_sh * sum c_sh(R2)  +  K4 * sum_{4a,4b} c_4(R2)/m^2  <=  theta * H,
H = f++ f+-/(f++ + f+-) (even J) or f+- (odd J);  theta = 1 (positivity) and theta = 1/2 (margin x2).
'b2' = nu^2 classes (B, C2) scaled by beta^2 (beyond-hybrid estimate, §3.4).  J* = smallest sampled J0 with every sampled J >= J0 passing.
Usage: python3 lsw_jstar.py -> cache/jstar_1000.txt
"""
import json, os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'cache', 'costs_1000.json')))
SC = {'A ': 86.75/6 - 8.0, 'B ': 11/6, 'C1': 1/6, 'C2': 1/6, 'C3': 11.75/6}
out = []; pr = lambda s: (print(s, flush=True), out.append(s))


def load(r, form, b2):
    c = r['cost']; b = r['beta']
    return sum(abs(SC[k])*abs(c[k][form])*(b*b if (b2 and k in ('B ', 'C2')) else 1) for k in SC), \
        sum(abs(c[k][1]) for k in c if k.startswith('sh')), sum(abs(c[k][1]) for k in ('4a', '4b'))/r['m']**2


def H(r):
    f = {h: v['floor'] for h, v in r['fl'].items()}
    return f['++']*f['+-']/(f['++'] + f['+-']) if len(f) == 2 else list(f.values())[0]


Js = sorted(set(r['J'] for r in R))
for form, b2, lab in ((0, False, 'R1 (de-shift), full CHMRS nu^2'), (0, True, 'R1, nu^2 classes x beta^2'),
                      (1, False, 'R2 (stripped), full'), (1, True, 'R2, nu^2 x beta^2')):
    for La, Ksh, K4 in ((1, 10, 1), (3, 10, 10), (10, 30, 10), (30, 100, 30), (100, 300, 100)):
        for th in (1.0, 0.5):
            fails = sorted(set(r['J'] for r in R if La*load(r, form, b2)[0] + Ksh*load(r, form, b2)[1] + K4*load(r, form, b2)[2] > th*H(r)))
            ok = [J for J in Js if all(J2 < J or J2 not in fails for J2 in Js)]
            pr('%-32s Lam_amp=%-4g K_sh=%-4g K4=%-4g theta=%.1f : J* = %-5s failing sampled J: %s' %
               (lab, La, Ksh, K4, th, (min(ok) if ok else 'none'), fails))
open(os.path.join(HERE, 'cache', 'jstar_1000.txt'), 'w').write('\n'.join(out) + '\n')
