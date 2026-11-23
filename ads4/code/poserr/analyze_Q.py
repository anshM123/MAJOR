"""
poserr/analyze_Q.py (POS_ERRORS.md §2.5): reads cache/q_large_m_<Dg>.json (q_large_m.py).  Per state (both channels):
   floor = 9/50 (CERT_CLOSE §2.3: certified for Q3* on J < 400, on b-hat <= 120, on ++ m <= 2M and on 2M <= m <= 10M, b-hat >= 120;
           float-checked >= 0.30 elsewhere) -- the flat FF/N itself is also printed
   cost  = 27/m (nu-flat order-1 factor on every image, sign-free) + g(J)*(Wq_lo + Wq_hi) + 0.6*Wd
           g(J) = 3 (450 + 100 J) c2 on the smeared k = 2, 3 parts (engine envelope E1', POS_ERRORS §2.3);
           |s_d| <= 0.6 (every measured channel 0.51 +- 0.08, m-independent from the gap to 20 Lam at J = 0, EESTEP §5)
   L     = 0.18/cost.  Also printed: Gq = 0.18/(Wq_lo + Wq_hi) (uniform sign-free tolerance on the smeared parts, c2) and the
           derivative-piece failure threshold s_fail' = 0.18/Wd.
Usage: python3 analyze_Q.py Dg
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); CA = os.path.join(HERE, 'cache')
BANDS = [0, 5, 10, 20, 40, 80, 160, 400.01]


def main(Dg):
    d = json.load(open(os.path.join(CA, 'q_large_m_%g.json' % Dg)))
    Lm = np.sqrt((Dg - 2)*(Dg - 1))
    rows = []
    for s in d['states']:
        J, mfl = s['J'], s['mfl']; m = mfl*Lm
        for h in ('++', '+-'):
            if h not in s:
                continue
            x = s[h]; g = 3*(450 + 100*J)
            D = (-(2*J + 1) + np.sqrt((2*J + 1)**2 + 4*m*m))/2 + J + 2
            cost = 27/m + g*(x['Wq_lo'] + x['Wq_hi']) + 0.6*x['Wd']
            rows.append(dict(J=J, mfl=mfl, bh=s['bhat'], h=h, FF=x['FF'], L=0.18/cost, Gq=0.18/(x['Wq_lo'] + x['Wq_hi']),
                             sfail=0.18/x['Wd'], fd=0.6*x['Wd']/cost, eng=(J*J/D <= 1 and J <= 141)))
    out = ['Q (Q3, delta = 0.2), Delta_g = %g: %d state-channels\n' % (Dg, len(rows))]
    for key, lab, fmt, agg in [('FF', 'min flat FF/N', '%8.4f', min), ('L', 'min leverage L = 0.18/cost', '%8.3g', min),
                               ('Gq', 'min Gq [c2]', '%8.2g', min), ('sfail', "min s_fail' = 0.18/Wd", '%8.3g', min),
                               ('fd', 'max share of cost from the derivative piece', '%8.3f', max)]:
        out.append('\n%s; rows m/Lam, columns b-hat bands %s\n' % (lab, ' '.join('[%g,%g)' % (BANDS[i], BANDS[i + 1]) for i in range(len(BANDS) - 1))))
        for mfl in sorted(set(r['mfl'] for r in rows)):
            c = []
            for i in range(len(BANDS) - 1):
                v = [r[key] for r in rows if r['mfl'] == mfl and BANDS[i] <= r['bh'] < BANDS[i + 1]]
                c.append(fmt % agg(v) if v else '%8s' % '-')
            out.append('  %5.1f  %s\n' % (mfl, ' '.join(c)))
    out.append('\nmin L on the engine-covered states (J^2/Delta <= 1, J <= 141) per m/Lam: ')
    for mfl in sorted(set(r['mfl'] for r in rows)):
        v = [r['L'] for r in rows if r['mfl'] == mfl and r['eng']]
        out.append('%g: %.3g (J <= %d);  ' % (mfl, min(v), max(r['J'] for r in rows if r['mfl'] == mfl and r['eng'])))
    out.append('\nmin Gq on the remaining states (J^2/Delta > 1 or J > 141) per m/Lam: ')
    for mfl in sorted(set(r['mfl'] for r in rows)):
        v = [r['Gq'] for r in rows if r['mfl'] == mfl and not r['eng']]
        out.append('%g: %.3g;  ' % (mfl, min(v)) if v else '%g: -;  ' % mfl)
    out.append('\n')
    for Lreq in (3.0, 10.0):
        out.append('\nL >= %g on every computed state with m/Lam >= ' % Lreq)
        for m0 in sorted(set(r['mfl'] for r in rows)):
            if all(r['L'] >= Lreq for r in rows if r['mfl'] >= m0):
                out.append('%g (all b-hat <= 400)' % m0); break
        else:
            out.append('none')
        out.append('\n')
    txt = ''.join(out)
    open(os.path.join(CA, 'analyze_Q_%g.txt' % Dg), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main(float(sys.argv[1]))
