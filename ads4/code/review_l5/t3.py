# Experiment 1: exact contour decomposition at infinite twist for the ideal psi2 weight.
import sys, numpy as np
from flint import arb, acb, ctx
from an_common import nf02, LAM2
ctx.prec = 128
LAM = arb(LAM2).sqrt(); Lf = float(LAM.mid())
p2, p3 = nf02()
P2 = [acb(arb(x.numerator)/x.denominator) for x in p2]
def psi(p):
    v = acb(0)
    for co in reversed(P2): v = v*p + co
    return v
I = acb(0,1)
def Fp(nu, beta):
    """outgoing (1-z)^{-3} P+_nu(cosh beta) = Gamma(i nu)/(sqrt(pi)Gamma(1/2+i nu)) e^{(-1/2+i nu)beta} 2F1(1/2,1/2-i nu;1-i nu;e^{-2beta})"""
    nu = acb(nu); beta = arb(beta)
    z = (beta/2).tanh()**2
    x = (-2*beta).exp()
    a = acb(0.5)-I*nu; cc = acb(1)-I*nu
    t = acb(1); h = acb(1); k = 0
    while True:
        t = t*(acb(0.5)+k)*(a+k)/((cc+k)*(k+1))*x; h += t; k += 1
        if k > 5 and abs(t).upper() < 1e-40*abs(h).lower(): break
    r = ((I*nu).lgamma() - (acb(0.5)+I*nu).lgamma()).exp()
    return r*((acb(-0.5)+I*nu)*beta).exp()*h/arb.pi().sqrt()*(1-z)**(-3)
def Freal(nu, beta):
    beta=arb(beta); z=(beta/2).tanh()**2
    return ((1-z)**(-3)*acb(beta.cosh()).legendre_p(acb(-0.5,nu),0,3)).real
def gl(f, a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    s = acb(0); sa = 0.0
    for xi, wi in zip(x, w):
        t = (a+b)/2 + (b-a)/2*xi
        v = f(t)*((b-a)/2*wi)
        s += v; sa += abs(complex(v.real.mid(), v.imag.mid()))
    return s, sa
beta = float(sys.argv[1]); delta = 0.25
# direct value G_A(c) = (2/Lam) int_0^Lam psi(nu/Lam) F(nu) dnu  (piecewise GL)
tot = 0.0; atot = 0.0
edges = np.linspace(0, Lf, 401)
for a_, b_ in zip(edges[:-1], edges[1:]):
    s, sa = gl(lambda t: acb(psi(arb(t)/LAM).real*2*Fp(t, beta).real), a_, b_, 24)
    tot += float(s.real.mid()); atot += sa
tot *= 2/Lf; atot *= 2/Lf
print(f'beta={beta}: G_A(c) direct = {tot:.10e}, int|.| = {atot:.3e}, d^3 G = {tot*float(((1-(arb(beta)/2).tanh()**2))**3):.4e}')
for Y in [2, 5, 10, 20, 30]:
    # bottom [0, delta]
    b0, ab0 = gl(lambda t: psi(acb(t)/LAM)*Fp(t, beta) if t > 0 else acb(0), 1e-9, delta, 40)
    # left vertical nu = delta + i y
    vl = acb(0); avl = 0.0
    ys = np.linspace(0, Y, int(Y*8)+1)
    for a_, b_ in zip(ys[:-1], ys[1:]):
        s, sa = gl(lambda y: psi((delta+I*y)/LAM)*Fp(delta+I*y, beta)*I, a_, b_, 20); vl += s; avl += sa
    # right vertical nu = Lam + i y
    vr = acb(0); avr = 0.0
    for a_, b_ in zip(ys[:-1], ys[1:]):
        s, sa = gl(lambda y: psi((LAM+I*y)/LAM)*Fp(LAM+I*y, beta)*I, a_, b_, 20); vr += s; avr += sa
    # top nu = t + iY, t in [delta, Lam]
    tp = acb(0); atp = 0.0
    ts = np.linspace(delta, Lf, 801)
    for a_, b_ in zip(ts[:-1], ts[1:]):
        s, sa = gl(lambda t: psi((t+I*Y)/LAM)*Fp(t+I*Y, beta), a_, b_, 16); tp += s; atp += sa
    k = 2/Lf
    tot2 = 2*(b0 + vl - vr + tp).real*k
    print(f' Y={Y}: sum = {float(tot2.mid()):.10e}; |pieces|: bottom {k*ab0:.2e} left {k*avl:.3e} right {k*avr:.2e} top {k*atp:.3e}; Re-values: left {float((2*k*vl.real).mid()):.4e} top {float((2*k*tp.real).mid()):.3e}')
