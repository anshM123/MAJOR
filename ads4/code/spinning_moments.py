r"""Diagonal spinning Regge kernels in d = 3: Chebyshev-Mehler formula and nu-moments (ads4/SPINNING_KERNELS.md).
  Phi^nu_{m,m}(b) = (sqrt2/pi) \int_0^b cos(nu t) T_{2m}(cosh(t/2)/cosh(b/2)) (cosh b - cosh t)^{-1/2} dt
  Q_0^{(m)} = T_{2m}(u0)/(2 sinh(b/2)),   u0 = sech(b/2)
  Q_2^{(m)} = -[ u0 T'_{2m}(u0)/(8 sinh(b/2)) + T_{2m}(u0)/(8 sinh^3(b/2)) ]
  Q_1^{(m)} = (sqrt2/pi) [ G(0)/b - \int_0^b (G(t)-G(0)) t^{-2} dt ],  G = T_{2m}(u(t)) (cosh b - cosh t)^{-1/2}
Also: derivative check  -(F'' + F'/sinh b)/(nu^2+1/4) = Phi_{1,1}  for F = P_{-1/2+i nu}(cosh b)  (a bi-vector component of
the scalar harmonic function is the weight-(1,1) matrix coefficient)."""
import mpmath as mp
mp.mp.dps = 20
T = lambda n, x: mp.chebyt(n, x)
dT = lambda n, x: n*mp.chebyu(n-1, x)
def Phi_diag(m, nu, b): return mp.re(mp.cosh(b/2)**(2*m)*mp.hyp2f1(0.5+m+1j*nu, 0.5+m-1j*nu, 1, -mp.sinh(b/2)**2))
def Phi_cheb(m, nu, b): return (mp.sqrt(2)/mp.pi)*mp.quad(lambda t: mp.cos(nu*t)*T(2*m, mp.cosh(t/2)/mp.cosh(b/2))/mp.sqrt(mp.cosh(b)-mp.cosh(t)), [0, b/2, b])
print("Chebyshev-Mehler formula vs hypergeometric form:")
for (m, nu, b) in ((4, 1.3, 0.8), (4, 5.0, 2.5), (2, 0.4, 3.0), (1, 2.2, 1.1)):
    print(f"  m={m} nu={nu} b={b}: {mp.nstr(Phi_cheb(m, mp.mpf(nu), mp.mpf(b)), 14)}   {mp.nstr(Phi_diag(m, mp.mpf(nu), mp.mpf(b)), 14)}")
F = lambda nu, b: mp.re(mp.hyp2f1(0.5+1j*nu, 0.5-1j*nu, 1, -mp.sinh(b/2)**2))
print("derivative check (weight (1,1)):")
for (nu, b) in ((1.3, 0.8), (3.0, 2.0)):
    nu, b = mp.mpf(nu), mp.mpf(b); F1 = mp.diff(lambda x: F(nu, x), b); F2 = mp.diff(lambda x: F(nu, x), b, 2)
    print(f"  nu={mp.nstr(nu,3)} b={mp.nstr(b,3)}: -(F''+F'/sinh)/(nu^2+1/4) = {mp.nstr(-(F2+F1/mp.sinh(b))/(nu**2+0.25), 12)}   Phi_11 = {mp.nstr(Phi_diag(1, nu, b), 12)}")
def Q0(m, b): return T(2*m, 1/mp.cosh(b/2))/(2*mp.sinh(b/2))
def Q2(m, b): u0 = 1/mp.cosh(b/2); return -(u0*dT(2*m, u0)/(8*mp.sinh(b/2)) + T(2*m, u0)/(8*mp.sinh(b/2)**3))
def Q1(m, b):
    def G(t):
        d = mp.cosh(b)-mp.cosh(t)
        return T(2*m, mp.cosh(t/2)/mp.cosh(b/2))/mp.sqrt(d) if d > 0 else mp.mpf(0)
    G0 = G(0)
    G2 = mp.diff(G, mp.mpf('1e-3'), 2)/2; t0 = b/100
    return (mp.sqrt(2)/mp.pi)*(G0/b - G2*t0 - mp.quad(lambda t: (G(t)-G0)/t**2, [t0, b/2, b]))
print("m=0 check: Q2(1) =", mp.nstr(Q2(0, mp.mpf(1)), 10), " -1/(8 sinh^3(1/2)) =", mp.nstr(-1/(8*mp.sinh(0.5)**3), 10))
print("small-beta limit of Q1: m=0:", mp.nstr(Q1(0, mp.mpf('0.05')), 6), "(1/24 = 0.041667);  m=4:", mp.nstr(Q1(4, mp.mpf('0.05')), 6), "(1/24 - 8 = -7.95833);  m=2:", mp.nstr(Q1(2, mp.mpf('0.05')), 6), "(1/24 - 2 = -1.95833)")
print("\nsigns of the moments in the helicity-flip-diagonal channel m = 4 (and m = 0 for comparison):")
print("  beta     Q0(m=0)   Q1(m=0)   Q2(m=0)  |   Q0(m=4)      Q1(m=4)      Q2(m=4)")
for b in map(mp.mpf, ('0.1', '0.2', '0.3', '0.4', '0.6', '0.8', '1.0', '1.5', '2.0', '3.0', '5.0')):
    print(f"  {mp.nstr(b,3):>4}  {mp.nstr(Q0(0,b),4):>9} {mp.nstr(Q1(0,b),4):>9} {mp.nstr(Q2(0,b),4):>9}  | {mp.nstr(Q0(4,b),5):>10} {mp.nstr(Q1(4,b),5):>12} {mp.nstr(Q2(4,b),5):>12}")
def zeros(f, lo, hi, n=4000):
    xs = [lo+(hi-lo)*mp.mpf(i)/n for i in range(n+1)]; z = []
    v = [f(x) for x in xs]
    for i in range(n):
        if v[i]*v[i+1] < 0: z.append(mp.findroot(f, (xs[i], xs[i+1]), solver='bisect'))
    return z
mp.mp.dps = 15
print("zeros in beta of Q0(m=4):", [mp.nstr(z, 6) for z in zeros(lambda b: Q0(4, b), mp.mpf('0.05'), mp.mpf(8), 800)])
print("zeros in beta of Q2(m=4):", [mp.nstr(z, 6) for z in zeros(lambda b: Q2(4, b), mp.mpf('0.05'), mp.mpf(8), 800)])
print("zeros in beta of Q1(m=4):", [mp.nstr(z, 6) for z in zeros(lambda b: Q1(4, b), mp.mpf('0.05'), mp.mpf(8), 120)])
