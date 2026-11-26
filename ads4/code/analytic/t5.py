# Experiment 2: finite-mass contour decomposition for the ideal psi2 weight; multiplier sizes on each piece.
import sys, time, numpy as np
from an_common import *
from kplus import *
ctx.prec = 3000
J,a,b = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); Ys = [float(x) for x in sys.argv[4].split(',')]
st = state(J,a,b); U = ladder(st)
beta = float((2*(st['c'].sqrt()).atanh()).mid()); d = st['d']
LAM = arb(LAM2).sqrt(); Lf = float(LAM.mid())
p2, p3 = nf02()
P2 = [acb(arb(x.numerator)/x.denominator) for x in p2]
def psi(p):
    v = acb(0)
    for co in reversed(P2): v = v*p + co
    return v
print(f'state J={J} a={a} b={b}: c={float(st["c"].mid()):.5f} beta={beta:.4f}')
def gl(f, a_, b_, n):
    x, w = np.polynomial.legendre.leggauss(n)
    s = acb(0); sa = 0.0; ra = []
    for xi, wi in zip(x, w):
        t = (a_+b_)/2 + (b_-a_)/2*xi
        v = f(t)
        s += v*((b_-a_)/2*wi); sa += abs(complex(v.real.mid(), v.imag.mid()))*((b_-a_)/2*wi)
    return s, sa
k = 2/Lf
stats = {}
def both(nu):
    Kf, _ = Kplus(nu, st, U); Ki, _ = Kplus(nu, st, U, inf=True)
    return Kf, Ki
t0 = time.time()
# real-axis exact value L(G_A^ideal) and infinite-twist value
def piece(name, fn, a_, b_, nseg, nn):
    edges = np.linspace(a_, b_, nseg+1)
    Sf = acb(0); Si = acb(0); Ad = 0.0; Ai = 0.0; mr = 0.0
    for e0, e1 in zip(edges[:-1], edges[1:]):
        x, w = np.polynomial.legendre.leggauss(nn)
        for xi, wi in zip(x, w):
            t = (e0+e1)/2 + (e1-e0)/2*xi; ww = (e1-e0)/2*wi
            nu, jac, wt = fn(t)
            Kf, Ki = both(nu)
            vf = wt*Kf*jac; vi = wt*Ki*jac
            Sf += vf*ww; Si += vi*ww
            Ad += float(abs(vf - vi).upper())*ww; Ai += float(abs(vi).upper())*ww
            r = float(abs(Kf - Ki).upper()/abs(Ki).lower()); mr = max(mr, r)
    print(f'  {name}: Re fm {float((2*k*Sf.real).mid()):+.6e}  Re inf {float((2*k*Si.real).mid()):+.6e}  2k*int|inf| {2*k*Ai:.3e}  2k*int|fm-inf| {2*k*Ad:.3e}  max|Kf/Ki-1| {mr:.3e}  [{time.time()-t0:.0f}s]', flush=True)
    return Sf, Si
delta = 0.25
tot_f = acb(0); tot_i = acb(0)
r = piece('bottom', lambda t: (acb(t), 1, psi(acb(t)/LAM)), 1e-12, delta, 1, 30); tot_f += r[0]; tot_i += r[1]
for Y in Ys:
    l = piece(f'left  Y={Y}', lambda y: (acb(delta, y), I, psi(acb(delta, y)/LAM)), 0, Y, int(4*Y)+1, 16)
    rr = piece(f'right Y={Y}', lambda y: (LAM + I*y, I, psi((LAM + I*y)/LAM)), 0, Y, int(2*Y)+1, 12)
    tp = piece(f'top   Y={Y}', lambda t: (acb(t, Y), 1, psi(acb(t, Y)/LAM)), delta, Lf, 400, 12)
    sf = tot_f + l[0] - rr[0] + tp[0]; si = tot_i + l[1] - rr[1] + tp[1]
    print(f' Y={Y}: total fm {float((2*k*sf.real).mid()):.8e}  total inf {float((2*k*si.real).mid()):.8e}')
