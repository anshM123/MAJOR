"""
h13_asym: large-impact-parameter asymptotics of the leading (Regge) D-safe action and check against the exact transforms.

Derivation [D] (h13 notes, NF02_LARGE_B.md section 3).  Mehler--Dirichlet: P_{-1/2+i nu}(cosh b) = (sqrt2/pi) int_0^b cos(nu t) K(t) dt,
K(t) = (cosh b - cosh t)^{-1/2}.  For phi(p) = sum_n c_n p^n (c_0 = 0 for NF02: psi(0) = 0) the endpoint (p = 0) expansion of
I = int_0^1 phi(p) P_{-1/2+i Lam p}(cosh b) dp at large Lam (fixed b) is
    I ~ sum_{n>=1} n! c_n M_n(b) / Lam^{n+1}  (+ p = 1 endpoint terms O(Lam^-6), phi ~ (1-p)^5),
    M_n = (sqrt2/pi) int_0^b K(t) [reg int_0^oo nu^n cos(nu t) dnu] dt :
      odd n : (sqrt2/pi) (-1)^{(n+1)/2} n! FP int_0^b K(t) t^{-n-1} dt      (non-local; vanishes in the flat limit),
      even n: (-1)^{n/2} K^{(n)}(0)/sqrt2                                   (local; flat limit = J0 Hankel moments).
  M_1 = (sqrt2/pi) g(b),  g(b) = K(0)/b - int_0^b (K(t)-K(0))/t^2 dt;  g -> sqrt2 pi/48 (b -> 0), so M_1 -> 1/24.
  M_2 = -(1/(2 sqrt2)) (cosh b - 1)^{-3/2}  (flat: -1/b^3).
For the D-safe action FR = sum_k u^{k-2} 2 eps_k sigma^k I_k:
    FR ~ A1 M_1/Lam^2 + A2 M_2/(2 Lam^3) + A3 M_3/(6 Lam^4) + ...,   A_n = sum_k u^{k-2} 2 eps_k sigma^k phi_k^{(n)}(0).
  NF02: A1 = 2 X(0) + sigma u Y(0) > 0, A2 = 2 psi2''(0) + sigma u psi3''(0) < 0 (psi3''(0) ~ 0 by construction), M_1 > 0 (checked here), M_2 < 0.
Output: cache/h13_asym.txt
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
from math import factorial
from h13_common import F, TR, KS, EPS, SIGMA, Lam_of, CACHE
mp.mp.dps = 40


def Kf(b, t):
    return (2*mp.sinh((b + t)/2)*mp.sinh((b - t)/2))**mp.mpf(-0.5)


def FPint(b, n):
    """Hadamard finite part of int_0^b K(t) t^{-n-1} dt, n odd (K even in t): subtract the Taylor polynomial of K to order n-1."""
    b = mp.mpf(b)
    ders = [mp.diff(lambda t: Kf(b, t), 0, j) for j in range(0, n + 2)]    # K^{(j)}(0), j <= n+1
    def integrand(s):
        t = b - s*s                                   # t = b - s^2 removes the (b-t)^{-1/2} endpoint singularity
        if t < mp.mpf(10)**-8*b:
            return 2*s*ders[n + 1]/mp.factorial(n + 1)
        K = (2*mp.sinh((b + t)/2)*mp.sinh(s*s/2))**mp.mpf(-0.5)
        tay = sum(ders[j]*t**j/mp.factorial(j) for j in range(0, n, 2))
        return 2*s*(K - tay)/t**(n + 1)
    I = mp.quad(integrand, [0, mp.sqrt(b)/4, mp.sqrt(b)])
    # finite part of int_0^b t^{j-n-1} dt for the subtracted Taylor terms (j - n - 1 != -1 since j even, n odd)
    for j in range(0, n, 2):
        I += ders[j]/mp.factorial(j)*b**(j - n)/(j - n)
    return I


def Mn(b, n):
    if n % 2:
        return float(mp.sqrt(2)/mp.pi*(-1)**((n + 1)//2)*mp.factorial(n)*FPint(b, n))
    return float((-1)**(n//2)*mp.diff(lambda t: Kf(mp.mpf(b), t), 0, n)/mp.sqrt(2))


def Acoef(h, u, n):
    s = SIGMA[h]
    return sum(u**(k - 2)*2*EPS[k]*s**k*factorial(n)*float(F.PHx[k][n] if n < len(F.PHx[k]) else 0) for k in KS)


if __name__ == '__main__':
    out = []
    def pr(s=''):
        print(s, flush=True); out.append(s)
    pr('# h13_asym: M_1 = (sqrt2/pi) g(beta) = int_0^oo 2 nu P_{-1/2+i nu}(cosh beta)/(e^{2 pi nu}+1) dnu (Mehler-Fock form, checked to 16 digits at beta = 0.3..14), M_2, M_3')
    for b in (1e-3, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 4.6, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0):
        M1, M2, M3 = Mn(b, 1), Mn(b, 2), Mn(b, 3)
        pr(f"  beta={b:7.3f}: M1 = {M1:+.6e} (flat-limit 1/24 = {1/24:.6e}), M2 = {M2:+.6e} (flat -1/b^3 = {-1/b**3:+.3e}), M3 = {M3:+.6e}")
    for b in (0.3, 1.0, 3.0, 6.0, 10.0, 14.0):
        mp.mp.dps = 25
        f = lambda nu: 2*nu*mp.re(mp.legenp(-0.5 + 1j*nu, 0, mp.cosh(b), type=3))/(mp.exp(2*mp.pi*nu) + 1)
        mf = float(mp.quad(f, [0, 0.5, 1, 2, 4, 8])); mp.mp.dps = 40
        pr(f"  Mehler-Fock check beta={b}: int 2nu P/(e^(2pi nu)+1) = {mf:.15e}  vs FP form {Mn(b, 1):.15e}")
    # elementary form: M_1 = (sqrt2/pi) int_0^b K(t) c(t) dt, c(t) = 2 int_0^oo nu cos(nu t)/(e^{2 pi nu}+1) dnu = cosh(t/2)/(4 sinh^2(t/2)) - 1/t^2
    def M1_elem(b):
        b = mp.mpf(b)
        c = lambda t: mp.cosh(t/2)/(4*mp.sinh(t/2)**2) - 1/t**2 if t > mp.mpf(10)**-6 else mp.mpf(1)/24
        f = lambda s: 2*s*c(b - s*s)*(2*mp.sinh((2*b - s*s)/2)*mp.sinh(s*s/2))**mp.mpf(-0.5)
        return float(mp.sqrt(2)/mp.pi*mp.quad(f, [0, mp.sqrt(b)/2, mp.sqrt(b)]))
    for b in (0.3, 3.0, 14.0, 20.0, 30.0, 40.0):
        v = M1_elem(b)
        pr(f"  elementary form beta={b}: M1 = {v:.15e}" + (f" (FP form {Mn(b, 1):.15e})" if b <= 14 else '') + f";  M1 e^(beta/2) beta/sqrt2 = {v*np.exp(b/2)*b/np.sqrt(2):.4f}")
    bs = np.concatenate([np.logspace(-3, 0, 40), np.linspace(1.0, 16.0, 151)])
    g = [Mn(b, 1) for b in bs]
    pr(f"  min M1 over {len(bs)} beta in [1e-3, 16]: {min(g):+.4e} at beta = {bs[int(np.argmin(g))]:.3f}  (positive: {min(g) > 0})")
    pr('')
    pr('# coefficients A_n = sum_k u^{k-2} 2 eps_k sigma^k phi_k^(n)(0)  (NF02 normalised, k(0) = 1)')
    for h in ('++', '+-'):
        for u in (1/2.25, 1/4, 1/9, 1/100, 0.0):
            pr(f"  {h} u={u:.4f}: A1 = {Acoef(h, u, 1):+.4f}  A2 = {Acoef(h, u, 2):+.4f}  A3 = {Acoef(h, u, 3):+.4f}  A4 = {Acoef(h, u, 4):+.4f}")
    pr('')
    pr('# FR (exact Mehler-Dirichlet transform) vs asymptotic A1 M1/L^2 + A2 M2/(2L^3) [+ A3 M3/(6 L^4)]')
    for Dg in (100, 300, 1000):
        L = Lam_of(Dg)
        for mfl in (1.5, 3.0, 10.0):
            u = 1/mfl**2
            for bh in (20, 50, 100, 200, 500, 1000):
                b = 2*np.arcsinh(bh/(2*L))
                if mfl*L*np.exp(-b/2) < Dg - 2:
                    continue
                M1, M2, M3 = Mn(b, 1), Mn(b, 2), Mn(b, 3)
                for h in ('++', '+-'):
                    I = TR.regge(L, b); s = SIGMA[h]
                    FR = sum(u**(k - 2)*2*EPS[k]*s**k*I[k] for k in KS)
                    a1 = Acoef(h, u, 1)*M1/L**2; a2 = Acoef(h, u, 2)*M2/(2*L**3); a3 = Acoef(h, u, 3)*M3/(6*L**4)
                    pr(f"  Dg={Dg:4d} m/L={mfl:4.1f} bhat~{bh:5d} beta={b:.4f} {h}: FR={FR:+.5e}  asym2={a1 + a2:+.5e} (rel {(a1 + a2)/FR - 1:+.1e})"
                       f"  asym3={a1 + a2 + a3:+.5e} (rel {(a1 + a2 + a3)/FR - 1:+.1e});  AdS/flat term ratio a1/a2 = {a1/a2:+.3f}")
    open(os.path.join(CACHE, 'h13_asym.txt'), 'w').write('\n'.join(out) + '\n')
