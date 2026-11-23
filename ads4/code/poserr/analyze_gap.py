"""
poserr/analyze_gap.py (POS_ERRORS.md §3): near-gap strip of Omega^EE at Delta_g = 1000 (cache/gap_m_1000.json from
`large_m_map.py 1000 gap 40`): per J, the minimum over the strip masses of the hybrid value, the sign-free tolerance Gtol (c2), the
leverage under the engine envelope E1' (analyze_D.budget) and the share of that cost carried by nu >= 20 (where, for J > 0, no
converged engine exists).  Also the calibration of the documented beyond-Regge corrections at m >= 1.3 Lam into effective uniform
|Gamma| (c2): Gamma_eff = |dF/N_s| / W  (EE_NEARGAP §6 values), using cache/large_m_<Dg>.json.
"""
import os, sys, json, glob
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); CA = os.path.join(HERE, 'cache')
sys.path.insert(0, HERE)
from analyze_D import budget


def main():
    d = json.load(open(os.path.join(CA, 'gap_m_1000.json'))); Lm = d['Lam']
    B = [budget(s, Lm) for s in d['states']]
    out = ['near-gap strip, Delta_g = 1000 (masses m_gap(J), 1.02, 1.05, 1.1, 1.2 Lam):\n',
           '  J | min lam_hyb [m/Lam] | min Gtol [c2] | min L_eng (E1\') | max share nu>=20 | J^2/Delta\n']
    for J in sorted(set(b['J'] for b in B)):
        bs = [b for b in B if b['J'] == J]
        i = int(np.argmin([b['lam'] for b in bs]))
        out.append('%3d | %8.4f [%6.4f] | %10.3g | %8.3g | %6.3f | %6.3f\n' % (
            J, bs[i]['lam'], bs[i]['mfl'], min(b['Gtol'] for b in bs), np.nanmin([b['L_eng'] for b in bs]),
            max(b['ftop'] for b in bs), max(b['j2d'] for b in bs)))
    # calibration: documented corrections -> effective uniform |Gamma|
    out.append('\ncalibration of documented beyond-Regge corrections (EE_NEARGAP §6) as uniform sign-free |Gamma| (c2):\n')
    cal = [(300, 1.3, J, v) for J, v in [(12, 0.04), (14, 0.07), (16, 0.12), (18, 0.19), (20, 0.26), (24, 0.46), (28, 0.49)]] + \
          [(300, 1.5, 24, 0.22), (300, 1.5, 28, 0.22), (300, 2.0, 24, 0.06), (300, 2.0, 28, 0.04), (1000, 1.3, 32, 0.005), (1000, 1.3, 40, 0.02)]
    for Dg, mfl, J, v in cal:
        f = os.path.join(CA, 'large_m_%g.json' % Dg)
        sts = [s for s in json.load(open(f))['states'] if abs(s['mfl'] - mfl) < 1e-9]
        if Dg == 1000:
            sts += [s for s in json.load(open(os.path.join(CA, 'large_m_1000b.json')))['states'] if abs(s['mfl'] - mfl) < 1e-9]
        s = min(sts, key=lambda s: abs(s['J'] - J))
        W = max(s[h]['W_lo'] + s[h]['W_hi'] for h in ('++', '+-') if h in s)
        g = 3*(450 + 100*s['J'])
        out.append('  Dg %4d, %.1f Lam, J = %2d (row J %2d): |dF/N| = %.3f -> Gamma_eff = %8.3g c2;  E1\' envelope %8.3g;  ratio %.2f\n' % (
            Dg, mfl, J, s['J'], v, v/W, g, v/W/g))
    txt = ''.join(out)
    open(os.path.join(CA, 'analyze_gap.txt'), 'w').write(txt)
    print(txt)


main()
