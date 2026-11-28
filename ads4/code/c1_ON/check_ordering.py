# Ordering test on S^1 x S^2 (beta=1): Y_{nJ} = (2w^2-u)/(w^2+u+1)^{5/2}, u = J(J+1)/r^2 (same conditionally
# convergent T_munu structure as delta).  D(r,Lam) = sum_{n,J}(2J+1) Y [exp(-P^2/Lam^2) - 1]  (covariant minus n-outer).
# Prediction: D -> -r^2 * 4pi * Delta_flat with Delta_flat = 1/(3 pi^2)... i.e. D/r^2 -> const, with NO r^0 term.
import numpy as np
def D(r, Lam):
    K = 1/r**2
    Jmax = int(40*Lam*r)
    J = np.arange(0, Jmax+1, dtype=float); u = J*(J+1)*K
    tot = 0.0
    nmax = int(8*Lam/(2*np.pi)) + 1
    for n in range(-nmax, nmax+1):
        w2 = (2*np.pi*n)**2
        Y = (2*w2 - u)/(w2 + u + 1)**2.5
        f = np.exp(-(w2 + u)/Lam**2)
        tot += ((2*J+1)*Y*(f - 1)).sum()
        # J-tail beyond Jmax of (2J+1)Y*(f-1) ~ -(2J+1)Y ~ +sum 2J u/u^{5/2}: add -int: (2J+1)(-u)/u^{5/2} -> -(2/K)... use integral
        uc = u[-1]; tot += -(-1)*(2/K)*(uc + w2 + 1)**-0.5 * 1.0  # int_{uc}^inf du/K * (-(-u)/(u+..)^{5/2}) approx
    return tot
for Lam in [20., 40., 80.]:
    vals = [D(r, Lam) for r in [2., 4., 8.]]
    rs = np.array([2., 4., 8.])
    # fit D = a r^2 + b
    A = np.vstack([rs**2, np.ones(3)]).T
    a, b = np.linalg.lstsq(A[:2], np.array(vals[:2]), rcond=None)[0]
    print(Lam, ['%.6f' % v for v in vals], 'a=%.6f b(r=2,4)=%.6f  pred a = -4/(3pi) * ... ' % (a, b), 'b from (4,8):',
          np.linalg.solve(A[1:], np.array(vals[1:]))[1])
print('-4/(3 pi) * (4 pi)/(4 pi^2)... flat: sum_n int d^2k/(2pi)^2 Y (f-1) -> -1/(3pi^2); times 4 pi r^2 =', -4*np.pi/(3*np.pi**2))
