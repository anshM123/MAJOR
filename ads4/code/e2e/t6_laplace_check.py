"""
Piecewise check of the tensor Fourier-Laplace transform (laplace.M_coeff) by direct numerical integration over H^2:
   Gamma(a) int d^2 yhat (-qhat.yhat)^{-a} (-yhat.l)^{-1/2-i nu} k_{q1q3}(yhat,l)[t]   vs   c(a, 1/2+i nu) (R + g1 S) e^{i mu theta0},
t = (z_1.eps_{h1})^l (z_3.eps_{h3})^l,  qhat = (1,0,0),  l = (1, cos theta0, sin theta0).
Includes entries with S != 0 (the parity-odd part, which is not probed by the parity-even partial waves of t4).
"""
import mpmath as mp, sympy as sp, itertools
from laplace import M_table, nu as symnu, a as syma
from formulas import gamma_a
mp.mp.dps = 18
out = []
def pr(*x):
    s = " ".join(str(v) for v in x); print(s, flush=True); out.append(s)
def md(u, v): return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def direct(l, q1, q3, h1, h3, nu, a, th0):
    lv = [mp.mpf(1), mp.cos(th0), mp.sin(th0)]
    eps = {1: [0, 1/mp.sqrt(2), 1j/mp.sqrt(2)], -1: [0, 1/mp.sqrt(2), -1j/mp.sqrt(2)]}
    s = mp.mpf(1)/2 + 1j*nu
    nth = 2*l + 1
    ths = [2*mp.pi*k/nth for k in range(nth)]
    def f(chi, phi):
        yh = [mp.cosh(chi), mp.sinh(chi)*mp.cos(phi), mp.sinh(chi)*mp.sin(phi)]
        yl = -md(yh, lv)
        n = [lv[i]/yl - yh[i] for i in range(3)]
        # n_perp_mu = eps_{mu a b} yh^a n^b ; contravariant: flip sign of component 0
        npc = [yh[1]*n[2] - yh[2]*n[1], yh[2]*n[0] - yh[0]*n[2], yh[0]*n[1] - yh[1]*n[0]]
        npv = [-npc[0], npc[1], npc[2]]
        def tz(h, q):
            vals = []
            for th in ths:
                z = [yh[i] + mp.cos(th)*n[i] + mp.sin(th)*npv[i] for i in range(3)]
                vals.append(md(z, eps[h])**l)
            return sum(vals[k]*mp.exp(-1j*q*ths[k]) for k in range(nth))/nth
        return mp.sinh(chi)*mp.cosh(chi)**(-a)*yl**(-s)*tz(h1, q1)*tz(h3, q3)
    val = mp.quad(f, [0, 0.7, 1.5, 3, 6, 14], mp.linspace(0, 2*mp.pi, 7))
    return mp.gamma(a)*val
pr("# direct H^2 integral  vs  c(a,s) (R + g1 S) e^{i mu theta0};   theta0 = 0.4")
th0 = mp.mpf('0.4')
tests = [(1, 1, -1, 1, -1, 1.3, 7), (1, 1, 0, 1, -1, 1.3, 7), (1, 1, 1, 1, 1, 0.6, 5), (1, 0, 1, -1, -1, 1.3, 7),
         (2, 2, -2, 1, -1, 1.3, 9), (2, 2, -1, 1, -1, 0.8, 9), (2, 1, 1, 1, 1, 1.3, 9), (2, 2, 2, -1, -1, 0.5, 11)]
import sys
if len(sys.argv) > 1: tests = tests[:int(sys.argv[1])]
for (l, q1, q3, h1, h3, nu, a) in tests:
    R, S, mu = M_table(l)[(q1, q3, h1, h3)]
    Rv = complex(sp.N(R.subs({symnu: sp.nsimplify(nu), syma: a}))) if R != 0 else 0
    Sv = complex(sp.N(S.subs({symnu: sp.nsimplify(nu), syma: a}))) if S != 0 else 0
    g1v = gamma_a(a - 1, nu)/gamma_a(a, nu)
    cas = 2**(a-1)*mp.sqrt(mp.pi)*gamma_a(a, nu)
    pred = cas*(Rv + g1v*Sv)*mp.exp(1j*mu*th0)
    dr = direct(l, q1, q3, h1, h3, mp.mpf(nu), a, th0)
    pr(f"  l={l} (q1,q3)=({q1},{q3}) (h1,h3)=({h1:+d},{h3:+d}) nu={nu} a={a}: direct = {mp.nstr(dr, 10)}   formula = {mp.nstr(pred, 10)}   [R = {Rv:.6g}, S = {Sv:.6g}, mu = {mu}]")
open('out_t6_laplace_check.txt', 'w').write("\n".join(out) + "\n")
