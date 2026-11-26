# exact real-axis value (2/Lam) int_0^Lam psi2(nu/Lam) K(nu) dnu with the terminating ladder [DV quadrature, exact integrand]
import sys, numpy as np
from an_common import *
ctx.prec = int(sys.argv[5]) if len(sys.argv) > 5 else 3000
J, a, b = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]); npan = int(sys.argv[4])
st = state(J, a, b); U = ladder(st)
LAM = arb(LAM2).sqrt(); Lf = float(LAM.mid())
p2, p3 = nf02(); P2 = [arb(x.numerator)/x.denominator for x in p2]
def psi(p):
    v = arb(0)
    for co in reversed(P2): v = v*p + co
    return v
x, w = np.polynomial.legendre.leggauss(16)
edges = np.linspace(0, Lf, npan + 1); tot = arb(0); tinf = arb(0)
c = st['c']; X = acb((1 + c)/(1 - c))
for e0, e1 in zip(edges[:-1], edges[1:]):
    for xi, wi in zip(x, w):
        t = (e0 + e1)/2 + (e1 - e0)/2*xi
        K = Kladder(3, arb(t), U)
        Ki = ((1 - c)**(-3)*X.legendre_p(acb(-0.5, t), 0, 3)).real
        ps = psi(arb(t)/LAM)
        tot += ps*K*((e1 - e0)/2*wi); tinf += ps*Ki*((e1 - e0)/2*wi)
print(f'real-axis: L(G) = {float((2/LAM*tot).mid()):.10e}   G(c) = {float((2/LAM*tinf).mid()):.10e}')
