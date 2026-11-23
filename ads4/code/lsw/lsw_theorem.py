"""
lsw/lsw_theorem.py (LARGESPIN_W.md §4-§5): evaluates the criterion of Theorem W on the sample states from cache/costs_1000.json.

For a state with per-helicity floors f_h and per-unit costs c_s (= |action of S/m^2|, s in a class set), uniform coefficient bound K:
  diagonal loss  d = K * sum_s c_s  (same for every entry, sign-free);
  2x2 (++, +-) Hermitian criterion, uniform in the coupling ratio (lambda_++, lambda_+-) in C^2:
      (f_++ - d)(f_+- - d) >= d^2   <=>   d <= f_++ f_+- / (f_++ + f_+-)
  1x1 (odd J, +-):  d <= f_+-.
  => K_tol(state, class set) = H / sum_s c_s,  H = f_++ f_+-/(f_++ + f_+-)  (even J)  or  f_+-  (odd J).
Class sets:
  S2   = {A, B, C1, C2, C3}                 O(1/m^2) CHMRS ring basis, de-shift form (R1)
  S2s  = same, stripped form (R2)           (worst case for R1)
  SH   = {sh+-1, sh+-2}                      O(1/m^2) index-shifted harmonics, stripped (they cannot be A-form)
  S4   = {4a, 4b} x 1/m^2                    O(1/m^4) representatives, stripped (worst)
Also: the scalar-CHMRS-shaped correction (exact (3.20) coefficients for Dphi = 3, k = 2 incl. the -8.0 P residual) and the
multiple of it tolerated; and the O(1/m) ring classes (for comparison: what parity buys).
Usage: python3 lsw_theorem.py -> cache/theorem_1000.txt
"""
import json, os, sys
sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'cache', 'costs_1000.json')))
out = []
pr = lambda s: (print(s, flush=True), out.append(s))

SC = {'A ': 86.75/6 - 8.0, 'B ': 11/6, 'C1': 1/6, 'C2': 1/6, 'C3': -11.75/6}      # scalar CHMRS (3.20) + residual, c2 units


def H(r):
    f = {h: v['floor'] for h, v in r['fl'].items()}
    if len(f) == 2:
        return f['++']*f['+-']/(f['++'] + f['+-']), f
    return list(f.values())[0], f


pr('J   m/L    beta  |  K_tol(S2,R1)  K_tol(S2,R2)  K_tol(SH)  K_tol(S4,R2) [c2^2-units: kappa_4/m^4]  | scalar-shape: used/H (R1, R2) |'
   ' O(1/m) ring A in (a/m) units: a_tol R1 / R2')
rows = []
for r in sorted(R, key=lambda q: (q['J'] % 2, q['J'], q['mfl'])):
    Hh, f = H(r)
    c = r['cost']; m = r['m']
    s2 = sum(abs(c[k][0]) for k in SC); s2s = sum(abs(c[k][1]) for k in SC)
    sh = sum(abs(c[k][1]) for k in c if k.startswith('sh'))
    s4 = sum(abs(c[k][1]) for k in ('4a', '4b'))/m**2
    scal1 = abs(sum(SC[k]*c[k][0] for k in SC)); scal2 = abs(sum(SC[k]*c[k][1] for k in SC))
    # O(1/m) ring amplitude (a/m) P with eta factor ~ 1: cost = a * m * c['A '] / cosh(beta)
    aA1 = Hh/(m*abs(c['A '][0])/np.cosh(r['beta'])); aA2 = Hh/(m*abs(c['A '][1])/np.cosh(r['beta']))
    row = dict(J=r['J'], mfl=r['mfl'], beta=r['beta'], K2=Hh/s2, K2s=Hh/s2s, KSH=Hh/sh, K4=Hh/s4, sc1=scal1/Hh, sc2=scal2/Hh,
               aA1=aA1, aA2=aA2)
    rows.append(row)
    pr('%3d %.3f %.4f | %10.3g %12.3g %10.3g %12.3g | %.1e  %.1e | %8.3g / %8.3g' %
       (r['J'], r['mfl'], r['beta'], row['K2'], row['K2s'], row['KSH'], row['K4'], row['sc1'], row['sc2'], aA1, aA2))

pr('')
pr('Minimum over sampled masses (m <= 1.3 Lam) for each J, and J*(K) = smallest sampled J0 such that every sampled state with J >= J0 '
   'has K_tol >= K (class sets as above):')
Js = sorted(set(q['J'] for q in rows))
mins = {J: {k: min(q[k] for q in rows if q['J'] == J) for k in ('K2', 'K2s', 'KSH', 'K4', 'aA1', 'aA2')} for J in Js}
for J in Js:
    pr('  J=%3d: K2 %.3g  K2s %.3g  KSH %.3g  K4 %.3g  aA(R1) %.3g  aA(R2) %.3g' % ((J,) + tuple(mins[J][k] for k in ('K2', 'K2s', 'KSH', 'K4', 'aA1', 'aA2'))))
even = [J for J in Js if J % 2 == 0]
for key, Ks in (('K2', [10, 30, 100, 300, 1000, 3000]), ('K2s', [10, 30, 100, 300, 1000]), ('KSH', [10, 100, 1000, 1e4]),
                ('K4', [1e3, 1e4, 1e5, 1e6])):
    for K in Ks:
        ok = [J for J in Js if all(mins[J2][key] >= K for J2 in Js if J2 >= J)]
        pr('  J*(%s, K=%g) = %s' % (key, K, min(ok) if ok else 'none sampled'))
open(os.path.join(HERE, 'cache', 'theorem_1000.txt'), 'w').write('\n'.join(out) + '\n')
json.dump(rows, open(os.path.join(HERE, 'cache', 'theorem_1000.json'), 'w'))
