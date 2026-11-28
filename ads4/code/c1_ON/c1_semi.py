"""Method B (semi-analytic): per-Matsubara-mode curvature coefficients for n != 0.
B_n = (1/3) delta0_n(0) + int_0^inf du D_n(u),
   delta0_n(u) = log(8 sqrt(w^2+u) Pi0_n(u)),  D_n(u) = -(u/3)[Pi0''/Pi0 - 3/(4 (w^2+u)^2)];
A_n = int_0^inf du delta0_n(u)   (flat r^2 coefficient, for the f_1 check).
"""
import sys, time
import flat
from flat import mp, mpf, log, sqrt, pi, quad, inf
mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 20

def Bn(n):
    w2 = (2*pi*n)**2
    d0 = log(8*sqrt(w2)*flat.Pi0(n, mpf(0)))
    def D(u):
        P = flat.Pi0(n, u); P2 = flat.Pi0(n, u, 2)
        return -(u/3)*(P2/P - mpf(3)/(4*(w2+u)**2))
    def dl(u):
        return log(8*sqrt(w2+u)*flat.Pi0(n, u))
    pts = [0, w2/4, w2, 4*w2, 16*w2, 100*w2, inf]
    I = quad(D, pts)
    A = quad(dl, pts)
    return d0, I, d0/3 + I, A

if __name__ == '__main__':
    n0, n1 = int(sys.argv[1]), int(sys.argv[2])
    for n in range(n0, n1+1):
        t = time.time()
        d0, I, B, A = Bn(n)
        print(n, mp.nstr(d0, 15), mp.nstr(I, 15), mp.nstr(B, 15), mp.nstr(A, 15), f'{time.time()-t:.1f}s', flush=True)
