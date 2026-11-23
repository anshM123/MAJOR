"""
n2_table.py (N2_SUBREGGE.md §2.3): corrected tolerance / leverage table at the H_LJ(1000) gap sample states, both floor conventions.
  * HLJ floor 0.9 lam0 (T_B; from cache/n2_tb.json, worse helicity +-)
  * LARGESPIN floor N_flat/5 (struct_1000.json, read-only): structure tolerances rescaled by exact grid cost ratios:
      G (window bump a e^{-(nu beta/3)^2}/m^2):  psi2 (LS)  |  psi2 + |X - psi2| (LS sec. 3.6)  |  |X| (sign-free, correct)
      A (ring (a/m) P(cosh beta)):  R1 = raw is A[env P] (de-shift identity: exact psi2 action = LS psi2 value)
                                    R2 = raw is env P itself: SIGNED exact action Re sum X env P  (vs LS's sign-free 0.31 at J = 344)
Usage: python3 n2_table.py -> cache/n2_table.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
import numpy as np
import n2_common as N

Dg = 1000; L = N.lam(Dg)
XD = N.xdata(Dg); nu, X, ps, dn = XD['nu'], XD['X'], XD['psi2'], XD['dnu']
ST = json.load(open(os.path.join(N.ADS4, 'code', 'largespin', 'cache', 'struct_1000.json')))
TB = {(r['J'], round(r['mfl'], 4)): r for r in json.load(open(os.path.join(N.CACHE, 'n2_tb.json')))}
lines = []


def pr(s=''):
    print(s, flush=True); lines.append(s)


pr('J  m/L   bhat | T_B(0.9lam0, c=3, +-): psi2 / |X| / |Y| | G a_tol (N_flat/5): psi2 / LS psi2+|X-psi2| / |X| | '
   'A a_tol (N_flat/5): R1(A-form, exact) / R2 signed exact / LS sign-free')
seen = set()
for q in ST:
    if q['reg'] != 'LJ' or q['J'] % 2 or q['J'] in seen:
        continue
    seen.add(q['J'])
    J, m, beta, mfl = q['J'], q['m'], q['beta'], q['mfl']
    h = q.get('h', '+-')
    envg = N.env('+-', nu/m, 2).real
    sel = nu < 0.999*m
    S = np.exp(-(nu*beta/3)**2)
    lo = nu < 30
    cps = (ps*S*envg).sum(); cX = (np.abs(X)*S*envg)[lo].sum() + (ps*S*envg)[~lo].sum(); cLS = cps + (np.abs(X - ps)*S*envg)[lo].sum()
    G = q['atol']['G']
    P = N.Pnu(beta, nu[sel]).real
    rA = (ps[sel]*P*envg[sel]).sum(); rA2 = ((X[sel]*P).real*envg[sel]).sum()
    A = q['atol']['A']
    tb = TB.get((J, round(mfl, 4)))
    tbs = ''
    if tb:
        t = tb['h']['+-']['3.0']
        tbs = f"{t['TB_psi']:6.0f} / {t['TB_X']:6.0f} / {t['TB_Y']:6.0f}"
    else:
        tbs = '     -  /     -  /     -'
    pr(f"{J:3d} {mfl:.3f} {q['bhat']:6.1f} | {tbs} | {G:6.0f} / {G*cps/cLS:6.0f} / {G*cps/cX:6.0f} | "
       f"{A:6.1f} / {A*abs(rA/rA2):7.3f} (R2/R1 = {rA2/rA:+8.2f}) / LS")
open(os.path.join(N.CACHE, 'n2_table.txt'), 'w').write('\n'.join(lines) + '\n')
