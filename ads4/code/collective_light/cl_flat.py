"""cl_flat: flat-space-limit model [M] of the collective O(1/c_T^2) light remainder of Omega^EE for <TTTT> (Einstein, AdS4).
Model: R = rho * |p_T omega_T|, with rho = one-loop/tree ratio of the smeared k = 2 arc at |s| ~ M^2 = Lam^2, 0 <= -t <= Lam^2:
   rho ~ c1 (G M^2/pi) log^2(Lam)     (D = 4 one-loop gravity: IR-divergent Weinberg term ~ (G s/pi) log(-s R^2)/eps, eps <-> 1/log Lam,
                                       plus finite (G/pi) s log^2 terms; AdS radius as IR regulator; c1 = O(1-10) not computed)
   G = 32/(pi c_T)   (supplement l.179: c_T = 32 L^2/(pi G), L = 1).
With X := pi c_T |R|/(512 n0 L_A) (m1_remainder) and |p_T omega_T| = 16 L_A n0 G one finds X = rho exactly, so the cubic bound needs rho <= 1.
This is an estimate of the size of a crossing-consistent one-loop light side; it is NOT a bound (README sec. 4).  Output: out/flat.txt
"""
import numpy as np, os
HERE = os.path.dirname(os.path.abspath(__file__))
L = ['# cl_flat: c_T* (rho = 1) for the flat-limit model rho = c1 (32 Lam^2/(pi^2 c_T)) log^2 Lam']
for Dg in (300, 1000):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    for c1 in (0.1, 1, 10):
        cts = c1*32*Lam**2/np.pi**2*np.log(Lam)**2
        L.append('  Delta_g=%5d  Lam=%7.1f  c1=%5.1f  ->  c_T* = %.1e   (log10 %.1f)' % (Dg, Lam, c1, cts, np.log10(cts)))
L.append('  compare: block-wise (m1_remainder) c_T* ~ 1e680-1e683 at Delta_g=1000; angular bound needs c_T >~ 2e8')
print('\n'.join(L)); open(os.path.join(HERE, 'out', 'flat.txt'), 'w').write('\n'.join(L) + '\n')
