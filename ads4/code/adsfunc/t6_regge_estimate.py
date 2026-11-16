"""t6: asymptotic Regge-regime estimate (beta >> 1/Lam) of the B3-vs-B2 competition in F_R (ADS4_FUNCTIONAL.md 3(b)):
   B2 ~ 2*65*Q1(beta)/Lam^2,  B3 ~ sigma u 162 Q0(beta)/Lam,  B4 = 10 W(beta) Lam^2 u^2;
   min_u (-162 Q0 u/Lam + 10 W Lam^2 u^2) = -(162 Q0)^2/(40 W Lam^4)  =>  B2 + that > 0  iff  Lam^2 > Lam_*^2 = (162 Q0)^2/(5200 W Q1).
Scalar moments (SPINNING_KERNELS.md Theorem K, m = 0): Q0 = 1/(2 sinh(beta/2)),
   Q1 = (sqrt2/pi)[G(0)/beta - int_0^beta (G(t)-G(0))/t^2 dt],  G(t) = (cosh beta - cosh t)^{-1/2}.  Floating point."""
import numpy as np
from scipy.integrate import quad
from common import curvT_all

def Q1(beta):
    G0 = (np.cosh(beta) - 1) ** -0.5
    G = lambda t: (np.cosh(beta) - np.cosh(t)) ** -0.5
    I = quad(lambda t: (G(t) - G0) / t ** 2, 0, beta, limit=400)[0]
    return np.sqrt(2) / np.pi * (G0 / beta - I)

print(' beta    Q0       Q1        W       Lam_*')
for beta in [0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]:
    Q0 = 1 / (2 * np.sinh(beta / 2))
    q1 = Q1(beta)
    W = curvT_all(1.0, beta)['W']
    print('%5.2f  %.4f  %.5f  %.5f  %7.2f' % (beta, Q0, q1, W, np.sqrt((162 * Q0) ** 2 / (5200 * W * q1))))
