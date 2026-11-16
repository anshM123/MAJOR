r"""Sharp threshold for positivity of the diagonal kernels Phi^nu_{m,m} in the complex spectral parameter nu.
(1) Orthogonality persists in the strip |Im nu| < j - 1/2:  \int Phi^nu_{m,m} cosh^{-2j}(b/2) sinh b db = 0, j <= m.
(2) For nu = i kappa with kappa >= m - 1/2, Phi^{i kappa}_{m,m}(b) > 0 for all b (Pfaff transformation: positive series).
(3) For nu = i kappa with kappa < m - 1/2 the kernel changes sign (it is orthogonal to a positive function)."""
import mpmath as mp
mp.mp.dps = 25
def Phi(m, nu, b): return mp.cosh(b/2)**(2*m)*mp.hyp2f1(0.5+m+1j*nu, 0.5+m-1j*nu, 1, -mp.sinh(b/2)**2)
print("(1) complex nu in the strip (m=4, j=4 needs |Im nu| < 3.5):")
for nu in (mp.mpc(1.0, 2.0), mp.mpc(0, 3.0), mp.mpc(2.5, 3.3), mp.mpc(0.7, 1.4)):
    I = mp.quad(lambda b: Phi(4, nu, b)*mp.cosh(b/2)**(-8)*mp.sinh(b), [0, 0.5, 1, 2, 4, 8, 16, 32, 64, 120, 250])
    A = mp.quad(lambda b: abs(Phi(4, nu, b))*mp.cosh(b/2)**(-8)*mp.sinh(b), [0, 0.5, 1, 2, 4, 8, 16, 32, 64, 120, 250])
    print(f"   nu = {nu}:  integral = {mp.nstr(I, 5)}   (integral of modulus {mp.nstr(A, 5)})")
print("(2),(3) sign of Phi^{i kappa}_{4,4}(b) on a grid b in (0, 12]:")
grid = [mp.mpf(i)/20 for i in range(1, 241)]
for kap in (0, 0.5, 1.5, 2.5, 3.0, 3.4, 3.5, 3.6, 4.5, 6.0):
    vals = [mp.re(Phi(4, 1j*mp.mpf(kap), b)) for b in grid]
    print(f"   kappa = {kap:>4}:  min = {mp.nstr(min(vals), 5):>12}   max = {mp.nstr(max(vals), 5):>12}   {'POSITIVE' if min(vals) > 0 else 'changes sign'}")
print("    kappa = 3.5: Phi = cosh^{-8}(b/2)?  ", mp.nstr(mp.re(Phi(4, 3.5j, mp.mpf(1.3))), 12), mp.nstr(mp.cosh(mp.mpf(1.3)/2)**(-8), 12))
print("(2') same threshold for m = 2 (kappa >= 1.5) and m = 0 (all kappa >= 0):")
for m, kaps in ((2, (1.0, 1.4, 1.5, 2.0)), (0, (0.0, 0.3, 1.0))):
    for kap in kaps:
        vals = [mp.re(Phi(m, 1j*mp.mpf(kap), b)) for b in grid]
        print(f"   m={m} kappa={kap}:  min = {mp.nstr(min(vals), 5)}   {'POSITIVE' if min(vals) > 0 else 'changes sign'}")
