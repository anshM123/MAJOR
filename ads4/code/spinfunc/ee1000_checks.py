"""
ee1000_checks (EE_DG1000_LARGEJ.md §1): (jacobi) scipy eval_jacobi vs the recurrence vs mpmath at complex c = 1 - 2(x + i/m)^2;
(deshift) the +- de-shift term with and without the fix on the (1040, 40) and (1040, 32) rows.
Usage: python3 ee1000_checks.py jacobi|deshift  -> cache/ee1000_checks_<mode>.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
mode = sys.argv[1]
out = []
def pr(s=''):
    print(s, flush=True); out.append(s)
if mode == 'jacobi':
    import mpmath as mp
    from scipy.special import eval_jacobi
    os.environ['EE1000_JACFIX'] = '1'
    import ee1000_common as E
    mp.mp.dps = 50
    pr("# |P^(0,8)_n(c)| at c = 1 - 2(x + i/m)^2, m = 1040: relative error of scipy eval_jacobi and of the recurrence vs mpmath (50 digits)")
    for n in (12, 20, 24, 28, 32, 36, 44, 52, 60):
        row = []
        for x in (0.3, 0.5, 0.6, 0.7, 0.8, 0.9):
            c = 1 - 2*(x + 1j/1040.0)**2
            ex = complex(mp.jacobi(n, 0, 8, mp.mpc(c.real, c.imag)))
            es = abs(eval_jacobi(n, 0, 8, c) - ex)/abs(ex); er = abs(E.jacobi_rec(n, 0, 8, c) - ex)/abs(ex)
            row.append(f"x={x}: scipy {es:.1e} rec {er:.1e}")
        pr(f"  n={n:2d}: " + "; ".join(row))
    for n in (40, 64):
        c = 1 - 2*(0.7 + 1j/1040.0)**2
        ex = complex(mp.legendre(n, mp.mpc(c.real, c.imag)))
        pr(f"  Legendre n={n}, x=0.7: numpy legval rel err {abs(np.polynomial.legendre.legval(c, np.eye(n+1)[n]) - ex)/abs(ex):.1e}, rec {abs(E.legendre_rec(n, c) - ex)/abs(ex):.1e}")
elif mode == 'deshift':
    os.environ['W2_DESPIKE'] = '0'
    import ee1000_common as E
    Q = E.Q
    for J, D in ((32, 1040.0), (40, 1040.0)):
        res = {}
        for fix in (False, True):
            if fix:
                Q.A.shift = Q._shift_generic(Q.A.env, E.Dc_fixed, Q.A.SIG)
            else:
                import l9_analysis as L9
                Q.A.shift = Q._shift_generic(Q.A.env, L9.Dc.__wrapped__ if hasattr(L9.Dc, '__wrapped__') else E._orig_Dc, Q.A.SIG)
            res[fix] = Q.A.grid(J, 2, ('C2', 'C3H'))[D]
        r0, r1 = res[False], res[True]
        pr(f"J={J} Delta={D:.0f}: a_h old {r0['a']} fixed {r1['a']}")
        for kd in ('C2', 'C3H'):
            for h in ('++', '+-'):
                pr(f"  {kd} {h} old  : " + " ".join(f"{v:9.3g}" for x, v in zip(r0['xs'], r0[(kd, h)]) if x >= 0.4))
                pr(f"  {kd} {h} fixed: " + " ".join(f"{v:9.3g}" for x, v in zip(r1['xs'], r1[(kd, h)]) if x >= 0.4))
open(os.path.join(HERE, 'cache', 'ee1000_checks_%s.txt' % mode), 'w').write("\n".join(out) + "\n")
