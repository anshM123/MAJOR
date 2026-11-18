"""Check h13_exact's rho against the cached e8 (Arb) records: same Delta (taken from the cache), mack_flint instead of mack_P_fast."""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import h13_exact as H
Dg, Df = int(sys.argv[1]), float(sys.argv[2])
db = pickle.load(open(os.path.join(HERE, 'cache', 'e1_%d_Df%g%s.pkl' % (Dg, Df, sys.argv[3] if len(sys.argv) > 3 else '')), 'rb'))
for key in [k for k in db if k[1] in (int(a) for a in sys.argv[4].split(','))][:4] if len(sys.argv) > 4 else list(db)[:2]:
    r0 = db[key]
    st = dict(J=key[1], mfl=r0['m']/H.Lam_of(Dg), Delta=r0['Delta'], eta=r0['eta'])
    r = H.rho_state(Dg, Df, st, 400)
    d = max(np.max(np.abs(r['rho'][k] - r0['rho'][k])) for k in H.KS)
    print(key, 'max |rho_h13 - rho_e8| =', d, ' arb rad', r['rad'], flush=True)
