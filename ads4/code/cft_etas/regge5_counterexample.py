"""Five-point Regge/superconvergence relations are LINEAR in the interference measure mu_n = g2(n) g3(n) (sign-indefinite).
Demo [P, linear algebra]: with heavy masses M_n >= Lambda, impose any finite set of relations sum_n mu_n M_n^{-2p} = L_p (p in P)
and show beta5 = sum_n mu_n M_n^{-2} (kernel independent of the constrained ones) takes any value, while the 4-point data
(g2(n), fixed) are untouched; g3(n) = mu_n/g2(n).  Contrast: a POSITIVE six-point moment bound sum g3^2/M^2 <= 2 C6 caps beta5 by
Cauchy-Schwarz, |beta5| <= 2 sqrt(A C6) with A = sum g2^2/(2 M^2)."""
import numpy as np
Lam = 1.0
M = Lam*np.array([1.0, 1.3, 1.7, 2.2, 3.0])
g2 = np.array([1.0, 0.8, 0.6, 0.5, 0.3])
P = [2, 3]                       # imposed five-point relations (say, set to zero = 'homogeneous superconvergence')
L = np.zeros(len(P))
for target in (-1e3, 0.0, 1e3):
    A = np.vstack([M**(-2*p) for p in P] + [M**-2.0])
    rhs = np.concatenate([L, [target]])
    mu, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    g3 = mu/g2
    print('beta5 target %+.0e -> achieved %+.3e ; relations residual %.1e ; six-point diagonal sum g3^2/M^2 = %.3e ; CS bound 2sqrt(A C6) = %.3e' % (
        target, np.sum(mu/M**2), np.max(np.abs(A[:-1] @ mu - L)), np.sum(g3**2/M**2),
        2*np.sqrt(np.sum(g2**2/(2*M**2))*np.sum(g3**2/(2*M**2)))))
