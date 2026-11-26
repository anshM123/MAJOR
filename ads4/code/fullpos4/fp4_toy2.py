"""
fullpos4/fp4_toy2.py (FULLPOS4.md §4) [DV]: end-to-end check of D + R (fp4_diag) against the exact ladder sum at a toy bandwidth with
strong cancellation.  G = 2 int psi2(p) F_{3, lam p} dp; exact L(G) = sum_k U_k M_k in mpmath (M_k by 200-node Gauss-Legendre in p, 60 digits).
Usage: fp4_toy2.py lam n N b xi0 [nA]
"""
import sys, os, math, json
sys.dont_write_bytecode = True
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fp4_int as I, fp4_diag as DG

lam = float(sys.argv[1]); n = int(sys.argv[2]); N = mp.mpf(sys.argv[3]); b = int(sys.argv[4]); xi0 = float(sys.argv[5])
nA = int(sys.argv[6]) if len(sys.argv) > 6 else 16
mp.mp.dps = 60
s = 2*b + 1; p = 3
w = I.Weight(2, lam)
xs, ws = mp.gauss_legendre if False else (None, None)
X, W = zip(*[(mp.mpf(x), mp.mpf(wt)) for x, wt in zip(*__import__('numpy').polynomial.legendre.leggauss(200))])
coef = [mp.mpf(c.numerator)/c.denominator for c in w.c]
M = [mp.mpf(0)]*(n + 1)
for x, wt in zip(X, W):
    pp = (x + 1)/2; ww = wt/2
    psi = mp.polyval(coef[::-1], pp)
    sn = (lam*pp)**2 + mp.mpf(1)/4
    g0, gm = mp.mpf(1), mp.mpf(0)
    for k in range(n + 1):
        M[k] += 2*psi*g0*ww
        g1 = ((2*k*k + 2*p*k + p - sn)*g0 - (k - 1 + p)**2*gm)/(k + 1)**2
        gm, g0 = g0, g1
U = mp.mpf(1); ex = mp.mpf(0); big = mp.mpf(0)
for k in range(n + 1):
    ex += U*M[k]; big = max(big, abs(U*M[k]))
    U *= (k - n)*(n + s + k)/((k - N)*(N + s + k))
a = float(N - n); J = n + b
a = int(a) if a == int(a) else a
r = DG.run(a, J, b, nA, lam, xi0)
out = dict(lam=lam, n=n, N=float(N), b=b, xi0=xi0, exact=float(ex), maxterm=float(big), D=r['D'], R=r['R_signed'], total=r['total'],
           relerr=abs(r['total'] - float(ex))/abs(float(ex)), Nrem_over=r['N_rem_soft']/abs(float(ex)), Dnorm_over=r['Dnorm']/abs(float(ex)))
print(json.dumps(out), flush=True)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'fp4_toy2.jsonl'), 'a').write(json.dumps(out) + '\n')
