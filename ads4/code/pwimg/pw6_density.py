"""
pwimg/pw6_density.py (PW_IMAGES.md §4.2): numerical check of the closed-form H2 densities of the index-shifted harmonics,
    Re P_{-1/2+i(nu-i)}(cosh b)  = cosh b P_nu(cosh b) - (1/2) int_0^b P_nu(cosh s) sinh s ds,
    Re P_{-1/2+i(nu-2i)}(cosh b) = cosh 2b P_nu(cosh b) - (1/2) int_0^b (3 cosh s + cosh b) P_nu(cosh s) sinh s ds,
(Abel inversion) at several (nu, b), mpmath.  Output: cache/pw6_density.txt
"""
import os, sys, io
sys.dont_write_bytecode = True
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 30
OUT = io.StringIO()
P = lambda mu, b: mp.legenp(-0.5 + 1j*mu, 0, mp.cosh(b), type=3)
worst = 0
for nu in (0.3, 1.7, 4.2, 9.0):
    for b in (0.4, 1.3, 3.0):
        l1 = (P(nu - 1j, b) + P(nu + 1j, b))/2
        r1 = mp.cosh(b)*P(nu, b) - 0.5*mp.quad(lambda s: P(nu, s)*mp.sinh(s), [0, b])
        l2 = (P(nu - 2j, b) + P(nu + 2j, b))/2
        r2 = mp.cosh(2*b)*P(nu, b) - 0.5*mp.quad(lambda s: (3*mp.cosh(s) + mp.cosh(b))*P(nu, s)*mp.sinh(s), [0, b])
        e1 = abs(l1 - r1)/abs(l1); e2 = abs(l2 - r2)/abs(l2); worst = max(worst, e1, e2)
        OUT.write(f'nu={nu} b={b}: shift 1 rel err {float(e1):.1e}, shift 2 rel err {float(e2):.1e}\n')
OUT.write(f'max rel err {float(worst):.1e}\n')
print(OUT.getvalue())
open(os.path.join(HERE, 'cache', 'pw6_density.txt'), 'w').write(__doc__ + '\n' + OUT.getvalue())
