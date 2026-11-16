r"""Positivity of the nu-moments Q_n(beta) of arXiv:2106.10274 (eqs. (5.43), (B.14), (B.15)) in d = 3.
(A) the paper's definition (B.15) at d = 3 (finite part):
      Q_n = -(sqrt2/pi) Gamma(n+1) sin(pi n/2) f.p. \int_0^beta t^{-n-1} (cosh beta - cosh t)^{-1/2} dt ;
(B) the formula of ads4/POSITIVITY_LEMMA.md:
      Q_n = (sqrt2/pi) \int_beta^infty K_n(t) (cosh t - cosh beta)^{-1/2} dt ,
      K_1(t) = -Im psi'(1 + i t/(2 pi)) / (2 pi^2)  > 0 ,   K_2(t) = -cosh(t/2)/(4 sinh^3(t/2)) < 0 ,   K_0 = coth(t/2)/2,
      K_n(t) = Gamma(n+1)[cos(pi n/2) t^{-n-1} + 2 sum_m sin((n+1) atan(t/(2 pi m))) ((2 pi m)^2+t^2)^{-(n+1)/2}]."""
import mpmath as mp
mp.mp.dps = 20
def fp_int(n, b):                         # f.p. \int_0^b t^{-n-1} g(t) dt for 0 < n < 2; g = g0 + g2 t^2 + O(t^4) near 0
    g = lambda t: (mp.cosh(b)-mp.cosh(t))**(-0.5); g0 = g(0); g2 = mp.mpf(1)/4*(mp.cosh(b)-1)**(-1.5)
    t0 = b/200                             # below t0 use the Taylor term (avoids cancellation in g(t) - g0)
    return g2*t0**(2-n)/(2-n) + mp.quad(lambda t: (g(t)-g0)*t**(-n-1), [t0, b/2, b]) - g0*b**(-n)/n
def QA(n, b):
    if n == 2: return -mp.sqrt(2)*mp.mpf(1)/4*(mp.cosh(b)-1)**(-1.5)       # limit n -> 2 of the definition: -1/(8 sinh^3(b/2))
    return -(mp.sqrt(2)/mp.pi)*mp.gamma(n+1)*mp.sin(mp.pi*n/2)*fp_int(n, b)
K = {0: lambda t: mp.coth(t/2)/2, 1: lambda t: -mp.im(mp.psi(1, 1+1j*t/(2*mp.pi)))/(2*mp.pi**2),
     2: lambda t: -mp.cosh(t/2)/(4*mp.sinh(t/2)**3)}
def Kgen(n, t, mmax=4000):
    s = mp.fsum(mp.sin((n+1)*mp.atan(t/(2*mp.pi*m)))*((2*mp.pi*m)**2+t**2)**(-(n+1)/2) for m in range(1, mmax))
    tail = (n+1)*t*(2*mp.pi)**(-n-2)*mmax**(-n-1)/(n+1)                         # leading tail estimate
    return mp.gamma(n+1)*(mp.cos(mp.pi*n/2)*t**(-n-1) + 2*(s+tail))
def QB(Kf, b): return (mp.sqrt(2)/mp.pi)*mp.quad(lambda t: Kf(t)/mp.sqrt(mp.cosh(t)-mp.cosh(b)), [b, b+1, b+5, b+20, b+80])
print("  n   beta      (A) definition         (B) formula            rel.diff")
for n in (1, 2):
    for b in map(mp.mpf, ('0.05', '0.5', '2', '6', '15')):
        a, c = QA(n, b), QB(K[n], b); print(f"  {n}  {mp.nstr(b,3):>5}  {mp.nstr(a, 13):>20} {mp.nstr(c, 13):>20}   {mp.nstr(abs(a-c)/abs(a), 3)}")
print("Q_0 from (B) vs 1/(2 sinh(beta/2)):", [(mp.nstr(QB(K[0], b), 12), mp.nstr(1/(2*mp.sinh(b/2)), 12)) for b in map(mp.mpf, ('0.5', '3'))])
print("Q_1(0.01) =", mp.nstr(QB(K[1], mp.mpf('0.01')), 10), "  1/24 =", mp.nstr(mp.mpf(1)/24, 10))
mp.mp.dps = 15
for n, b in ((mp.mpf('0.5'), mp.mpf(1)), (mp.mpf('0.5'), mp.mpf(4))):
    a = QA(n, b); c = (mp.sqrt(2)/mp.pi)*mp.quad(lambda t: Kgen(n, t, 600)/mp.sqrt(mp.cosh(t)-mp.cosh(b)), [b, b+1, b+6, b+40], maxdegree=5)
    print(f"  n=0.5 beta={mp.nstr(b,2)}:  (A) {mp.nstr(a,9)}   (B, series) {mp.nstr(c,9)}   rel.diff {mp.nstr(abs(a-c)/abs(a),2)}")
print("K_1 closed form vs series at t=1.7:", mp.nstr(K[1](mp.mpf('1.7')), 10), mp.nstr(Kgen(mp.mpf(1), mp.mpf('1.7'), 3000), 10))
print("min K_1 on t in (0, 60]:", mp.nstr(min(K[1](mp.mpf(t)/10) for t in range(1, 601)), 5), "  max K_2:", mp.nstr(max(K[2](mp.mpf(t)/10) for t in range(1, 601)), 5))
