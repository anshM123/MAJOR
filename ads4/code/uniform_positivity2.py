r"""G_M[f](beta) = \int_0^1 f(p) P_{-1/2+iMp}(cosh beta) dp = (sqrt2/pi) \int_0^beta ftil(M t) (cosh beta - cosh t)^{-1/2} dt,
ftil(tau) = \int_0^1 f(p) cos(p tau) dp (closed form via sympy).  Grid scan in (M, beta) for the graviton weight f_g3 and
comparison with (i) the flat-space functional  fhat(b) = \int_0^1 f(p) J_0(p b) dp  at b = M beta, and
(ii) the two-term curved-regime expansion a_1 Q_1(beta)/M^2 + a_2 Q_2(beta)/M^3."""
import sympy as sp, mpmath as mp
mp.mp.dps = 18
p, tau = sp.symbols('p tau', positive=True)
f_sym = (1-p)**3*p*(65+155*p+47*p**2)
ft_sym = sp.simplify(sp.integrate(sp.expand(f_sym)*sp.cos(p*tau), (p, 0, 1)))
ser = sp.series(ft_sym, tau, 0, 8).removeO()
ft_exact = sp.lambdify(tau, ft_sym, 'mpmath'); ft_ser = sp.lambdify(tau, ser, 'mpmath')
ftil = lambda x: ft_ser(x) if x < 0.5 else ft_exact(x)
f = sp.lambdify(p, f_sym, 'mpmath')
def G(M, b):
    n = max(2, int(M*b/3)+1); pts = [b*mp.mpf(i)/n for i in range(n)] + [b*(1-mp.mpf(1)/(4*n)), b]
    def integrand(t):
        d = mp.cosh(b)-mp.cosh(t)
        return ftil(M*t)/mp.sqrt(d) if d > 0 else mp.mpf(0)
    return (mp.sqrt(2)/mp.pi)*mp.quad(integrand, pts)
fhat = lambda x: mp.quad(lambda q: f(q)*mp.besselj(0, q*x), mp.linspace(0, 1, int(x/3)+2))
K1 = lambda t: -mp.im(mp.psi(1, 1+1j*t/(2*mp.pi)))/(2*mp.pi**2)
Q1 = lambda b: (mp.sqrt(2)/mp.pi)*mp.quad(lambda t: K1(t)/mp.sqrt(mp.cosh(t)-mp.cosh(b)), [b, b+1, b+5, b+20, b+80])
Q2 = lambda b: -1/(8*mp.sinh(b/2)**3)
print("check G against direct definition (M=5, beta=0.7):", mp.nstr(G(5, mp.mpf('0.7')), 10),
      mp.nstr(mp.quad(lambda q: f(q)*mp.re(mp.hyp2f1(0.5+5j*q, 0.5-5j*q, 1, -mp.sinh(mp.mpf('0.35'))**2)), mp.linspace(0, 1, 6)), 10))
a1, a2 = 65, -40
for M in (5, 20, 80):
    print(f"M = {M}:    beta      M*beta      G_M(beta)      flat fhat(M beta)       [a1 Q1/M^2 + a2 Q2/M^3]")
    worst = None
    for x in (0.5, 1, 2, 3, 5, 8, 12, 20, 40, 80, 160, 400):
        b = mp.mpf(x)/M
        if b > 12: continue
        g = G(M, b); fl = fhat(mp.mpf(x)); cur = a1*Q1(b)/M**2 + a2*Q2(b)/M**3
        print(f"          {mp.nstr(b,4):>8} {x:>8}   {mp.nstr(g, 6):>14}   {mp.nstr(fl, 6):>14}        {mp.nstr(cur, 6):>14}")
    grid = [mp.mpf(k)/(2*M) for k in range(1, 12*M if M <= 20 else 240)]
    vals = [(G(M, b), b) for b in grid]; mn = min(vals, key=lambda v: v[0])
    print(f"   min over beta grid (step {mp.nstr(mp.mpf(1)/(2*M),3)}, up to {mp.nstr(grid[-1],3)}): G = {mp.nstr(mn[0], 6)} at beta = {mp.nstr(mn[1], 4)}")
