"""
fullpos4/fp4_toy.py (FULLPOS4.md §2.2) [DV]: end-to-end check of the frequency-integrated outgoing §55 contour at a small state
with a scaled bandwidth Lam_toy:  G(z) = 2 int_0^1 psi2(p) F_{3, Lam_toy p}(z) dp,  exact L(G) = sum_{k<=n} U_k M_k,
M_k = 2 int psi2(p) g_k(Lam_toy p) dp (Gauss-Legendre, 400 nodes).   Usage: fp4_toy.py lam n N b
"""
import sys, os, math
sys.dont_write_bytecode = True
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fp4_int as I, fp4_abs as AB

lam = float(sys.argv[1]); n = int(sys.argv[2]); N = float(sys.argv[3]); b = int(sys.argv[4])
s = 2*b + 1; p = 3
w = I.Weight(2, lam)
x, wx = np.polynomial.legendre.leggauss(400); pp = (x + 1)/2; wp = wx/2
psi = w.psi(pp); sn = (lam*pp)**2 + 0.25
g = [np.ones_like(pp)]; gm = np.zeros_like(pp)
for k in range(n + 1):
    g.append(((2*k*k + 2*p*k + p - sn)*g[k] - (k - 1 + p)**2*(g[k - 1] if k else 0))/(k + 1)**2)
M = [2*(psi*g[k]*wp).sum() for k in range(n + 1)]
U = 1.0; ex = 0.0; big = 0.0
for k in range(n + 1):
    ex += U*M[k]; big = max(big, abs(U*M[k]))
    U *= (k - n)*(n + s + k)/((k - N)*(N + s + k))
a = N - n; J = n + b
r = AB.run(a, J, b, amp=I.Amp(2, lam=lam), toy=True)
print(f'lam={lam} n={n} N={N} b={b}: exact L(G)={ex:.12e} (max|U_k M_k|={big:.2e})  contour={r["signed"]:.12e}  '
      f'rel.err={abs(r["signed"]-ex)/abs(ex):.2e}  int|W|/|L(G)|={r["Nabs"]/abs(ex):.2e}', flush=True)
