"""repro_engine/dimdiff.py -- convergence diagnostics: Gamma_K per harmonic (c2) of direct runs relative to a base run, and to the engine.
Usage: python3 dimdiff.py Delta base_ns ns1 ns2 ...   (ns as 8-8-8-8-6-8-8; f64 runs; numax from the files)
"""
import sys, os
import numpy as np
import compare as CM
import sdchain as C

SHOWNU = (0.5, 0.9, 1.3, 1.5, 2.02, 2.5, 3.02, 3.5, 4.02, 6.02, 8.02, 10.02)


def gk(D, ns, prec='f64'):
    nus, rho = CM.direct_ratio(D, ns, prec)
    g = C.gammaK(nus, rho, D)
    return nus, g['GK'].real/g['sabs']


if __name__ == '__main__':
    D = float(sys.argv[1]); base = sys.argv[2]; others = sys.argv[3:]
    prec = os.environ.get('PREC', 'f64')
    en, er = C.engine_ratio(D, 5); ge = C.gammaK(en, er, D); e5 = ge['GK'].real/ge['sabs']
    nb, gb = gk(D, base, prec)
    idx = [int(round((v - 0.02)/0.04)) for v in SHOWNU if v < nb[-1]]
    print('nu:            ' + ' '.join(f'{nb[i]:8.2f}' for i in idx))
    print(f'base-K5 {base}: ' + ' '.join(f'{gb[i] - e5[i]:+8.3f}' for i in idx))
    for ns in others:
        prec2 = prec
        if ':' in ns:
            ns, prec2 = ns.split(':')
        n2, g2 = gk(D, ns, prec2)
        k = min(len(n2), len(nb))
        print(f'{ns}:{prec2} - base: ' + ' '.join(f'{g2[i] - gb[i]:+8.3f}' for i in idx if i < k) +
              '   | -K5: ' + ' '.join(f'{g2[i] - e5[i]:+8.3f}' for i in idx if i < k))
