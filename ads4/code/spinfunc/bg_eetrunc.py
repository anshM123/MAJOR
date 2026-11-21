"""
bg_eetrunc (BELOW_FOLD2.md §3.4): truncation table of the EE images: per-order sizes m^2 [log R_{<=k} - log R_{<=k-1}] (relative to
nu_0 = 0.02), R = I_st/I_ref, for EE k = 4 and k = 2 (state order k, reference at K_ref = 2).
Usage: python3 bg_eetrunc.py Delta J K -> appends cache/bg_eetrunc.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import bg_eeq as EQ

NUP = (0.5, 1.5, 2.02, 2.5, 3.02, 4.02, 5.02, 8.02, 10.02, 13.02, 18.02, 25.02, 30.02, 40.02)

if __name__ == '__main__':
    D, J, K = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    db = pickle.load(open(EQ.fn(D, J, K), 'rb'))
    m = np.sqrt((D - J - 2)*(D + J - 1))
    nus = db['st']['nus']
    idx = [int(np.argmin(np.abs(nus - v))) for v in NUP]
    lines = [f'## Delta={D} J={J} m={m:.2f}: per-order m^2[log R_k - log R_(k-1)] - (nu_0), k = 1..{K}, at nu = ' + ' '.join(f'{v:g}' for v in NUP)]
    for kind in ('EE4', 'EE2'):
        Iref, _ = EQ.value(db['ref'], kind)
        prev = None
        for Kt in range(0, K + 1):
            Ist, _ = EQ.value(db['st'], kind, Kt=Kt)
            L = m*m*np.log(Ist/Iref)
            if prev is not None:
                d = L - prev; d = d - d[0]
                lines.append(f'   {kind} k={Kt}: ' + ' '.join(f'{d[i].real:+.2f}' for i in idx))
            prev = L
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'bg_eetrunc.txt'), 'a').write('\n'.join(lines) + '\n')
