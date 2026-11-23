"""
largespin/ls_strip.py (LARGESPIN_ANALYTIC.md §5): H_bulk as a relative statement, and the strip on which a stated bulk factor suffices.

Reads (read-only) POS_ERRORS' per-state weights code/poserr/cache/large_m_<Dg>*.json (W_lo: nu < 12, W_mid: 12-20, W_top: nu >= 20,
off-diagonal analogues; flat N, computed curvature, Lemma-C floor 'rig').  Envelope E1' exactly as poserr/analyze_D.budget:
nu-flat 27/m (x1), non-flat g(J) = 3(450 + 100 J) c2 on nu < 20 (engine, [D given RD] + NE, x3).  For nu >= 20 (the bulk, J > 0) the
envelope is f * g(J) and we print, per state, the largest f for which the state stays positive (the "bulk factor"):
    f_bulk = (floor - 27/m - g (W_lo + W_mid)) / (g W_top)        (per helicity, worst; off-diagonal included)
floor = rig*N_s (certified N/5 + Lemma C, used where b-hat < 10) or N_f/5 + min(0, curv) (computed curvature, b-hat >= 10).
In relative terms (correction / flat envelope E_k(x)) the bulk tolerance is  eps_bulk = f_bulk * g(J) / m^2.
Output: cache/strip_<Dg>.txt
"""
import os, sys, json, glob
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
PE = os.path.join(HERE, '..', 'poserr', 'cache')
Dg = float(sys.argv[1]) if len(sys.argv) > 1 else 1000.0


def delta_of(m, J):
    al = (-(2*J + 1) + np.sqrt((2*J + 1)**2 + 4*m*m))/2
    return al + J + 2


rows = []
for f in sorted(glob.glob(os.path.join(PE, 'large_m_%g*.json' % Dg))):
    d = json.load(open(f)); Lm = d['Lam']
    for s in d['states']:
        J, mfl = s['J'], s['mfl']; m = mfl*Lm; D = delta_of(m, J)
        g = 3*(450 + 100*J); eps = 27/m
        best = None
        for h in ('++', '+-'):
            if h not in s:
                continue
            q = s[h]
            fl_rig = q['rig']                               # already / N_s
            fl_cmp = (q['Nf']/5 + min(0.0, q['cv']))/q['Ns']
            fl = fl_rig if s['bhat'] < 10 else fl_cmp
            lo = g*(q['W_lo'] + q['W_mid'] + s['W_off_lo'] + s['W_off_mid'])
            top = g*(q['W_top'] + s['W_off_top'])
            fb = (fl - eps - lo)/top if top > 0 else np.inf
            L1 = fl/(eps + lo + top)
            r = dict(h=h, fl=fl, fb=fb, L1=L1, share=top/(eps + lo + top))
            if best is None or fb < best['fb']:
                best = r
        rows.append(dict(J=J, mfl=mfl, bh=s['bhat'], j2d=J*J/D, m=m, g=g, **best, eps_rel=best['fb']*g/m**2))
B = rows
out = []
w = out.append
w('Delta_g = %g: %d states; floor = Lemma-C rig (b < 10) / computed (b >= 10); E1\' with bulk factor f on nu >= 20' % (Dg, len(B)))
bands = [0, 5, 7.5, 10, 15, 20, 30, 40]
w('\nmin bulk factor f_bulk (states with J^2/Delta <= 1 only = engine-covered below nu = 20); rows m/Lam, cols b-hat bands %s' % bands)
for mfl in sorted(set(b['mfl'] for b in B)):
    row = []
    for i in range(len(bands) - 1):
        v = [b['fb'] for b in B if b['mfl'] == mfl and bands[i] <= b['bh'] < bands[i + 1] and b['j2d'] <= 1]
        row.append('%9.3g' % min(v) if v else '%9s' % '-')
    w('  %5.1f %s' % (mfl, ' '.join(row)))
w('\nsame, relative bulk tolerance eps_bulk = f_bulk g(J)/m^2 (fraction of the flat envelope E_k(x) the nu >= 20 correction may reach)')
for mfl in sorted(set(b['mfl'] for b in B)):
    row = []
    for i in range(len(bands) - 1):
        v = [b['eps_rel'] for b in B if b['mfl'] == mfl and bands[i] <= b['bh'] < bands[i + 1] and b['j2d'] <= 1]
        row.append('%9.3g' % min(v) if v else '%9s' % '-')
    w('  %5.1f %s' % (mfl, ' '.join(row)))
for freq in (76, 30, 10, 3):
    w('\nlargest b-hat per m/Lam such that EVERY state below it (J^2/Delta <= 1) has f_bulk >= %g:' % freq)
    line = []
    for mfl in sorted(set(b['mfl'] for b in B)):
        bs = sorted((b['bh'], b['fb'], b['j2d']) for b in B if b['mfl'] == mfl)
        bm = 0.0
        for bh, fb, jd in bs:
            if jd > 1 or not (fb >= freq):
                break
            bm = bh
        Jm = int((bm*mfl - 1)/2)
        line.append('%g: b < %.1f (J <= %d)' % (mfl, bm, Jm))
    w('  ' + ';  '.join(line))
txt = '\n'.join(out)
print(txt)
open(os.path.join(HERE, 'cache', 'strip_%g.txt' % Dg), 'w').write(txt + '\n')
