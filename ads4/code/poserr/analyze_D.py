"""
poserr/analyze_D.py (POS_ERRORS.md §2.2-2.3): reads cache/large_m_<Dg>*.json (large_m_map.py) and prints, per m/Lam and b-hat band,
  lam_hyb  : min over helicity of the hybrid value (exact functions, quadrature)                                   [NE]
  rig      : min over helicity of (N_flat/5 - |curvature| bound)/N_s  (certified NF02 floor + Lemma C)             [certified + P]
  Gtol     : uniform sign-free |Gamma| (c2, every nu, diagonal and off-diagonal) that can make the state fail, using rig
  L_eng    : rig / cost under the engine envelope E1 (below), only where J^2/Delta <= 1 (engine convergent regime)
  L_HR     : rig / cost under the explicit sign-free hypothesis H_R': |Gamma| <= 30 m + 1e4 c2 on every nu, nu-flat 3*27/m
Envelope E1 (POS_ERRORS §2.3, from eng_orders.py; x3 safety on every term):
  nu-flat (common) part of order 1: relative 27/m (EE4; EE2 21/m) on the whole image -> cost <= 3*27/m in lam/N_s (sign-free);
  non-flat, nu < 12:  g_lo(J)   = 3 (450 + 100 J) c2      (Delta-independent at fixed J: J = 0 at Delta = 300, 1e3, 2e4; J = 28..100)
  non-flat, nu >= 12: g_hi(J,m) = 3 (450 + 0.3 J m) c2    (grows ~ m at fixed J: J = 28..100 at Delta = 1100, 2e4)
  k = 3 (C3H) part: the same envelope [M] (no C3H engine); its share of the cost is printed.
Usage: python3 analyze_D.py Dg
"""
import os, sys, json, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); CA = os.path.join(HERE, 'cache')
BANDS = [0, 5, 10, 20, 40, 80, 160, 400.01]


def delta_of(m, J):
    al = (-(2*J + 1) + np.sqrt((2*J + 1)**2 + 4*m*m))/2
    return al + J + 2


def budget(s, Lm):
    J, mfl = s['J'], s['mfl']; m = mfl*Lm; D = delta_of(m, J)
    hs = [h for h in ('++', '+-') if h in s]
    lam = min(s[h]['lam'] for h in hs); rig = min(s[h]['rig'] for h in hs)
    # floor2: certified flat floor N_f/5 plus the COMPUTED (exact-function) curvature where it is negative  [certified + NE]
    fl2 = {h: (s[h]['Nf']/5 + min(0.0, s[h]['cv']))/s[h]['Ns'] for h in hs}
    floor2 = min(fl2.values())
    Wtot = max(s[h]['W_lo'] + s[h]['W_hi'] for h in hs) + s['W_off_lo'] + s['W_off_hi']
    Gtol = floor2/Wtot if floor2 > 0 else -1.0
    # engine envelope E1' (x3 safety); windows nu < 12 | 12-20 | >= 20 (the last is NOT engine-covered for J > 0)
    eps = 27/m; g = 3*(450 + 100*J)          # nu-flat factor: measured to 1e-3 (x1); non-flat engine terms x3
    cost, sh_top, sh_k3 = 0.0, 0.0, 0.0
    for h in hs:
        lo = g*(s[h]['W_lo'] + s['W_off_lo']); mid = g*(s[h]['W_mid'] + s['W_off_mid']); top = g*(s[h]['W_top'] + s['W_off_top'])
        c = eps + lo + mid + top
        if c/fl2[h] > cost/max(floor2, 1e-300) or cost == 0.0:
            cost = c; sh_top = top/c; sh_k3 = g*(s[h]['W3_lo'] + s[h]['W3_hi'] + s['W3_off'])/c
    L_eng = min(fl2[h] for h in hs)/max(eps + g*(s[h]['W_lo'] + s[h]['W_hi'] + s['W_off_lo'] + s['W_off_hi']) for h in hs)
    gg = 30*m + 1e4
    L_hr = floor2/max(eps + gg*(s[h]['W_lo'] + s[h]['W_hi'] + s['W_off_lo'] + s['W_off_hi']) for h in hs)
    return dict(J=J, mfl=mfl, bh=s['bhat'], D=D, j2d=J*J/D, lam=lam, rig=rig, floor2=floor2, Gtol=Gtol,
                L_eng=(L_eng if (J*J/D <= 1 and J <= 141) else np.nan), f3=sh_k3, ftop=sh_top, L_HR=L_hr,
                cvrel=max(s[h]['cvb']/s[h]['Nf'] for h in hs))


def main(Dg):
    sts = []
    for f in glob.glob(os.path.join(CA, 'large_m_%g*.json' % Dg)):
        d = json.load(open(f)); Lm = d['Lam']; sts += d['states']
    B = [budget(s, Lm) for s in sts]
    out = []
    w = out.append
    w('Delta_g = %g, Lam = %.4f, %d states\n' % (Dg, Lm, len(B)))
    for key, lab, fmt in [('lam', 'min lam_hyb', '%8.4f'), ('rig', 'min rig floor (cert. flat/5 - Lemma C)', '%8.4f'),
                          ('floor2', 'min floor2 (cert. flat/5 + computed curvature if negative)', '%8.4f'),
                          ('Gtol', 'min Gtol [c2] (uniform sign-free |Gamma| on every nu that can make the state fail; from floor2)', '%8.2g'),
                          ('L_eng', 'min leverage, engine envelope E1 (J^2/Delta <= 1 only; nan = no such state)', '%8.3g'),
                          ('f3', 'max share of the E1 cost from the k = 3 part', '%8.3f'),
                          ('ftop', 'max share of the E1 cost from nu >= 20 (not engine-covered for J > 0)', '%8.3f'),
                          ('L_HR', "min leverage under H_R' (|Gamma| <= 30m + 1e4, all nu)", '%8.3g')]:
        w('\n%s; rows m/Lam, columns b-hat bands %s\n' % (lab, ' '.join('[%g,%g)' % (BANDS[i], BANDS[i + 1]) for i in range(len(BANDS) - 1))))
        for mfl in sorted(set(b['mfl'] for b in B)):
            row = []
            for i in range(len(BANDS) - 1):
                v = [b[key] for b in B if b['mfl'] == mfl and BANDS[i] <= b['bh'] < BANDS[i + 1] and np.isfinite(b[key])]
                row.append((fmt % (max(v) if key in ('f3', 'ftop') else min(v))) if v else '%8s' % '-')
            w('  %5.1f  %s\n' % (mfl, ' '.join(row)))
    # largest b-hat (per m) up to which every state has L >= Lreq
    for key in ('L_eng', 'L_HR'):
        for Lreq in (3.0, 10.0):
            w('\n%s >= %g holds for every state with b-hat <= bmax(m/Lam):  ' % (key, Lreq))
            for mfl in sorted(set(b['mfl'] for b in B)):
                bs = sorted((b['bh'], b[key]) for b in B if b['mfl'] == mfl)
                bm = 0.0
                for bh, L in bs:
                    if not (np.isfinite(L) and L >= Lreq):
                        break
                    bm = bh
                w('%g: %.1f;  ' % (mfl, bm))
            w('\n')
    txt = ''.join(out)
    open(os.path.join(CA, 'analyze_D_%g.txt' % Dg), 'w').write(txt)
    print(txt)
    json.dump(B, open(os.path.join(CA, 'analyze_D_%g.json' % Dg), 'w'))


if __name__ == '__main__':
    main(float(sys.argv[1]))
