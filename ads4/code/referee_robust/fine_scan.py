"""Referee: fine u-scan (m in [1,40]) of FF/N, J <= 600, to locate true minima (floating point)."""
import sys, numpy as np
from heavy_indep import weights
def run(name, nq=700):
    psi2, psi3, a = weights(name)
    P2 = np.polynomial.polynomial.Polynomial([float(c) for c in psi2]); P3 = np.polynomial.polynomial.Polynomial([float(c) for c in psi3])
    xg, wg = np.polynomial.legendre.leggauss(nq); pp = (xg + 1) / 2; ww = wg / 2
    ms = np.unique(np.concatenate([np.linspace(1, 6, 1500), np.linspace(6, 40, 6000)]))
    for h in ('++', '+-'):
        sig = 1 if h == '++' else -1; c4 = 1 if h == '++' else 41; best = (9, None)
        for m in ms:
            u = 1 / m**2; x = 1 - 2*u*pp*pp; nmax = 600
            D = np.empty((nmax+1, nq)); D[0] = 1
            if h == '++':
                D[1] = x
                for n in range(2, nmax+1): D[n] = ((2*n-1)*x*D[n-1] - (n-1)*D[n-2]) / n
                Js = np.arange(nmax+1)
            else:
                D[1] = 1 + 10*(x-1)/2
                for n in range(2, nmax+1):
                    c = 2*n + 8.0
                    D[n] = ((c-1)*(c*(c-2)*x - 64)*D[n-1] - 2*(n-1)*(n+7)*c*D[n-2]) / (2*n*(n+8)*(c-2))
                Js = np.arange(nmax+1) + 4
            B2 = D @ (ww*P2(pp)*(2-u*pp*pp)); B3 = sig*u*(D @ (ww*P3(pp))); FW = float(a)*(2*Js*(Js+1)-c4+sig)*u**3
            if h == '++':
                s = Js % 2 == 0; Js, B2, B3, FW = Js[s], B2[s], B3[s], FW[s]
            r = (B2+B3+FW)/(abs(B2)+abs(B3)+abs(FW)); i = np.argmin(r)
            if r[i] < best[0]: best = (r[i], (round(m, 4), int(Js[i]), round((2*Js[i]+1)/m, 3)))
        print(name, h, 'fine min FF/N = %.5f at (m, J, bhat) = %s' % best, flush=True)
for n in sys.argv[1:]: run(n)
