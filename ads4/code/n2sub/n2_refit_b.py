"""n2_refit_b.py (N2_SUBREGGE.md §3.4): F1 of n2_refit.py plus the CHMRS ring-derivative class nu^2 D (B; raw A-form), at the graviton
window rows: is the non-ring remainder r compatible with zero once the full ring class (A: D, B: nu^2 D, C: x D') is allowed?"""
import os, sys, pickle
sys.dont_write_bytecode = True
import numpy as np
import n2_common as N
from n2_refit import Acol, dprime, lsq, d
lines = []
for key in [(100, 1100.0, 12.0, 34.0, 'cv', ('_w',)), (100, 1100.0, 0.0, 34.0, 'cv', ('_none',)), (64, 1064.0, 12.0, 34.0, 'cv', ('_w',)),
            (64, 364.0, 0.0, 34.0, 'cv', ('_w',))]:
    r = d[key]; J = key[0]; m = r['m']; nus = np.asarray(r['nus']); sel = (nus >= 12) & (nus <= 34.0001); nu = nus[sel]; x = nu/m
    for h in ('++', '+-'):
        y = np.asarray(r[('C2', h)]['y'])[sel]; err = np.asarray(r[('C2', h)]['err'])[sel]
        AD = Acol(lambda z: N.Dsh(h, J, z), nu, m, h); ADp = Acol(lambda z: dprime(h, J, z)*z/(2*J + 1), nu, m, h)
        AB = Acol(lambda z: (z*m)**2*N.Dsh(h, J, z)/100.0, nu, m, h); A1 = Acol(lambda z: 0*z + 1.0, nu, m, h)
        dk = Acol(lambda z: N.kappa('C2', h, J, m, z), nu, m, h) - N.kappa('C2', h, J, m, x).real
        if len(nu) >= 5:
            c, e, rms, res = lsq([AD, AB, ADp, A1], y - dk)
            chi = np.sqrt(np.mean((res/np.maximum(err, 1))**2))
            s = (f"{key[1]:g},{J} {key[5][0]} {h}: c0 {c[0]:+.0f}  b(nu^2/100) {c[1]:+.0f} +- {e[1]:.0f}  s {c[2]:+.0f}  r {c[3]:+.0f} +- {e[3]:.0f}  "
                 f"rms {rms:.0f} (rms/err {chi:.2f}); nodes {len(nu)}")
        else:
            s = f"{key}: too few nodes"
        print(s); lines.append(s)
open(os.path.join(N.CACHE, 'n2_refit_b.txt'), 'w').write('\n'.join(lines) + '\n')
