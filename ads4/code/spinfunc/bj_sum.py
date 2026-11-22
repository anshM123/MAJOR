"""bj_sum: per-state lam/N_s (J, m/Lam) from hc_chain z4 pickles.  Usage: python3 bj_sum.py J tag1 tag2 ..."""
import os, sys, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
J = int(sys.argv[1])
rows = {}
for t in sys.argv[2:]:
    d = pickle.load(open(os.path.join(HERE, 'cache', 'hc_chain_z4_%s.pkl' % t), 'rb'))
    for r in d[300]:
        if r['J'] == J:
            rows.setdefault(r['mfl'], {})[t] = float(r['lam'])
print('   m/Lam  ' + ' '.join(f'{t[-22:]:>22s}' for t in sys.argv[2:]))
for mf in sorted(rows):
    print(f'   {mf:5.3f}  ' + ' '.join(f'{rows[mf].get(t, np.nan):+22.4f}' for t in sys.argv[2:]))
print('   min    ' + ' '.join(f'{min(rows[m].get(t, 9) for m in rows):+22.4f}' for t in sys.argv[2:]))
