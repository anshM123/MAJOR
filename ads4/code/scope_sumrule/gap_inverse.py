"""
Inverse (gap) statement from the thermal sum rules.  [computed]

tau_*(D) = 160/(9(D-2)) [1 - 49/(2D+9)^2]^{-1/2}   (letter.tex eq. (angular)).
Even channel alone:  delta_E = |S_E - 3/8| = |t4|/32 <= tau/32.
With the Hall channel (if |S_H| = |alpha_T|/32):  delta_th = sqrt((S_E-3/8)^2 + S_H^2) = tau/32.
tau <= tau_*(D)  <=>  (D-2) sqrt(1 - 49/(2D+9)^2) <= 5/(9 delta).
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

tau_star = lambda D: 160/(9*(D - 2))/np.sqrt(1 - 49/(2*D + 9)**2)
lhs = lambda D: (D - 2)*np.sqrt(1 - 49/(2*D + 9)**2)

# monotonicity of tau_* on D > 2
Ds = np.linspace(2.0001, 1e5, 200001)
print('tau_* strictly decreasing on (2, 1e5]:', bool(np.all(np.diff(tau_star(Ds)) < 0)))
# c_max = sup_{D>2} [ (D-2) - lhs(D) ]  =>  D <= 2 + c_max + 5/(9 delta)
res = minimize_scalar(lambda D: -((D - 2) - lhs(D)), bounds=(2.0, 50.0), method='bounded')
cmax = -res.fun
print('sup_D [(D-2) - (D-2)sqrt(1-49/(2D+9)^2)] = %.6f at D = %.4f' % (cmax, res.x))
print('=> D_g <= %.4f + 5/(9 delta)  (collaborator: 2.26 + 5/(9 delta))' % (2 + cmax))
# the collaborator's 2.26: check it is a valid (slightly loose) constant
print('2 + c_max <= 2.26:', 2 + cmax <= 2.26)


def g_star_exact(delta):
    """largest D with tau_*(D) >= 32 delta (tau_* decreasing => D_g <= this)."""
    t = 32*delta
    return brentq(lambda D: tau_star(D) - t, 2 + 1e-12, 1e12)


print('\n delta      exact D_max     2 + 5/(9 delta)   2.26 + 5/(9 delta)')
for dl in (1/8, 0.1, 0.03, 0.01, 1e-3, 1e-4, 1e-6):
    print('%8.2e  %14.6f  %16.6f  %16.6f' % (dl, g_star_exact(dl), 2 + 5/(9*dl), 2.26 + 5/(9*dl)))
print('\nNote: |t4| <= 4 (collider) gives delta_E <= 1/8, and tau <= 4 gives delta_th <= 1/8.')
print('At Delta_g = 1000: tau_*/32 = %.4e (max allowed deviation of S_E from 3/8 without positivity);' % (tau_star(1000)/32))
print('                   2.89e-8/32 = %.2e (with positivity, no light single traces).' % (2.89e-8/32))
