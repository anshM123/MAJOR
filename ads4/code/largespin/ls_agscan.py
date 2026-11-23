"""
largespin/ls_agscan.py (LARGESPIN_ANALYTIC.md §4.3): the travelling-wave (AG) saddle test of HLJ_ANALYTIC §3, applied to H_R' states
(m >= 1.3 Lam, large b-hat) instead of the near-gap H_LJ states.  Thin wrapper: imports spinfunc/lj_scan.py and lj_saddle.py read-only
(their functions write nothing) and writes only to largespin/cache/agscan_<Dg>.txt.
Per state: min_nu sgn*Re kappa~/beta on nu >= 2m/(2J+1) for both waves (local frequency = decay rate in Im nu), the smallest base
Hessian eigenvalue (non-degeneracy), and the decay rate / beta along vertical lines nu_r + i y, y <= 30/beta.
Usage: python3 ls_agscan.py Dg m/Lam:J,m/Lam:J,...
"""
import os, sys
os.environ.setdefault('L1_DEG', '2')
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
import numpy as np
import lj_scan as LS

Dg = float(sys.argv[1])
specs = [(float(a), int(b)) for a, b in (s.split(':') for s in sys.argv[2].split(','))]
Lam = LS.Lam_of(Dg)
f = open(os.path.join(HERE, 'cache', 'agscan_%d.txt' % Dg), 'a')


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); f.write(s + '\n'); f.flush()


for mL, J in specs:
    mgap = np.sqrt((Dg - 2.0)*(Dg - 1.0 + 2*J))
    mfac = mL*Lam/mgap
    r = LS.scan_state(Dg, J, mfac, pr=pr)
    pr('Dg=%g J=%d m/Lam=%.4f (mfac=%.4f) Delta=%.3f m=%.2f beta=%.5f bhat=%.2f' % (Dg, J, r['m']/Lam, mfac, r['Delta'], r['m'], r['beta'],
                                                                                 (2*J + 1)*Lam/r['m']))
    for sgn in (1, -1):
        rows = r['res'][sgn]
        okall = all(x[1] for x in rows)
        nus = np.array([x[0] for x in rows]); kr = np.array([x[2] for x in rows])
        sel = nus >= 2*r['m']/(2*J + 1)
        kmin = np.min(sgn*kr[sel]) if sel.any() else np.nan
        numin = nus[sel][np.argmin(sgn*kr[sel])] if sel.any() else np.nan
        eb = np.min([x[5] for x in rows]); nub = nus[np.argmin([x[5] for x in rows])]
        pr('  wave %+d: Newton ok=%s nu in [%.1f, %.1f]; min sgn Re kappa~/beta (nu >= %.1f) = %.5f at nu=%.1f; min base |eig| = %.3g at nu=%.1f'
           % (sgn, okall, nus.min(), nus.max(), 2*r['m']/(2*J + 1), kmin/r['beta'], numin, eb, nub))
    nuB = 5*r['m']/(2*J + 1)
    for nu_r in (0.9*Lam, 0.5*Lam, 0.1*Lam, nuB):
        try:
            beta, out = LS.vertical(Dg, J, mfac, nu_r, [1/r['beta'], 3/r['beta'], 10/r['beta'], 30/r['beta']])
            pr('  vertical nu_r=%.1f: ' % nu_r + '  '.join('y=%.2f ok=%d rate/beta=%.4f' % o for o in out))
        except Exception as e:
            pr('  vertical nu_r=%.1f: failed (%s)' % (nu_r, e))
f.close()
