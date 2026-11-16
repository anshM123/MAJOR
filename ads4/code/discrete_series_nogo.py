r"""Orthogonality of the principal-series kernels Phi^nu_{m,m} to the lowest-weight discrete-series coefficient
cosh^{-2m}(beta/2)  (Haar measure sinh(beta) d beta on K\G/K-type functions of weight (m,m)):
     I_m(nu) := \int_0^infty Phi^nu_{m,m}(beta) cosh^{-2m}(beta/2) sinh(beta) d beta = 0   for every real nu, m >= 1.
Consequence: A(beta) = \int F(nu) Phi^nu_{m,m}(beta) d nu cannot be >= 0 for all beta unless it vanishes identically.
For m = 0 there is no discrete series and Phi^nu_{0,0} = P_{-1/2+i nu} has no such constraint."""
import mpmath as mp
mp.mp.dps = 25
def Phi(m, nu, b): return mp.re(mp.cosh(b/2)**(2*m)*mp.hyp2f1(0.5+m+1j*nu, 0.5+m-1j*nu, 1, -mp.sinh(b/2)**2))
print("  m    nu      I_m(nu)            integral of |integrand| (for scale)")
for m in (1, 2, 4):
    for nu in map(mp.mpf, ('0', '0.4', '1.3', '3.7')):
        f = lambda b: Phi(m, nu, b)*mp.cosh(b/2)**(-2*m)*mp.sinh(b)
        pts = [0, 0.5, 1, 2, 4, 8, 16, 30] if m > 1 else [0, 0.5, 1, 2, 4, 8, 16, 30, 60, 90]
        I = mp.quad(f, pts); A = mp.quad(lambda b: abs(f(b)), pts)
        print(f"  {m}   {mp.nstr(nu,3):>4}   {mp.nstr(I, 6):>14}        {mp.nstr(A, 6)}")
print("check with the lowest K-type removed (should NOT vanish): \\int Phi^nu_{4,4} cosh^{-6}(b/2) sinh b db, nu=1.3:",
      mp.nstr(mp.quad(lambda b: Phi(4, mp.mpf('1.3'), b)*mp.cosh(b/2)**(-6)*mp.sinh(b), [0, 0.5, 1, 2, 4, 8, 16, 30]), 6))
print("m = 0 (no constraint): \\int P_{-1/2+i nu} cosh^{-8}(b/2) sinh b db, nu=1.3:",
      mp.nstr(mp.quad(lambda b: Phi(0, mp.mpf('1.3'), b)*mp.cosh(b/2)**(-8)*mp.sinh(b), [0, 0.5, 1, 2, 4, 8, 16, 30]), 6))
