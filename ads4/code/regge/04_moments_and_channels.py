r"""Check 4: (a) flat-space benchmark: the stripped helicity partial waves of arXiv:2201.06602 (A.3),
    dtilde^J_{4,4}(x) = 2F1(4-J, J+5; 1; (1-x)/2),  dtilde^J_{0,0}(x) = P_J(x),  x = 1 - 2p^2/m^2,
coincide at fixed b = 2J/m up to O(1/m^2) (both -> J_0(bp), their (3.5), (A.11));
(b) moments Q_n^{(m)}(beta) = \int_0^infty nu^n Phi^nu_{m,m}(beta) dnu (Abel-regularised), n = 0,1,2, m = 0,2,4:
closed forms of SPINNING_KERNELS.md Theorem K against an independent direct nu-integration of the 2F1 form;
(c) the quantities relevant for the corrected channel assignment (REGGE_FROM_BLOCKS.md, section 6)."""
import mpmath as mp
mp.mp.dps = 20
print("(a) flat benchmark, b = 3, p = 0.7:   m, dtilde_44 - J_0(bp), dtilde_00 - J_0(bp), m^2 (dtilde_44 - dtilde_00)")
b, p = mp.mpf(3), mp.mpf('0.7')
for m in (20, 40, 80, 160):
    J = int(b*m/2); z = p**2/mp.mpf(m)**2
    d44 = mp.hyp2f1(4-J, J+5, 1, z); d00 = mp.hyp2f1(-J, J+1, 1, z); j0 = mp.besselj(0, b*p)
    print(f"   m={m:4d} J={J:4d}: {mp.nstr(d44-j0, 6):>12} {mp.nstr(d00-j0, 6):>12} {mp.nstr(m**2*(d44-d00), 8):>14}")
# ---------------- (b) moments ----------------
def T2m(m, u): return mp.chebyt(2*m, u)
def dT2m(m, u): return mp.diff(lambda x: mp.chebyt(2*m, x), u)
def Q0(m, be):
    u0 = 1/mp.cosh(be/2); return T2m(m, u0)/(2*mp.sinh(be/2))
def Q2(m, be):
    u0 = 1/mp.cosh(be/2); sh = mp.sinh(be/2)
    return -(u0*dT2m(m, u0)/(8*sh) + T2m(m, u0)/(8*sh**3))
def Q1(m, be):
    """(sqrt2/pi) [G(0)/beta - int_0^beta (G(t)-G(0)) t^-2 dt], G(t) = T_2m(cosh(t/2)/cosh(beta/2)) (cosh beta - cosh t)^(-1/2);
    substitution t = beta cos(tau) removes the endpoint singularity."""
    def G_stable(tau):                   # returns G(t) with cosh(be) - cosh(t) computed without cancellation
        t = be*mp.cos(tau); dt = 2*be*mp.sin(tau/2)**2
        return T2m(m, mp.cosh(t/2)/mp.cosh(be/2))/mp.sqrt(2*mp.sinh((be + t)/2)*mp.sinh(dt/2)), t
    G0 = T2m(m, 1/mp.cosh(be/2))/mp.sqrt(mp.cosh(be) - 1)
    G2 = mp.diff(lambda t: T2m(m, mp.cosh(t/2)/mp.cosh(be/2))/mp.sqrt(mp.cosh(be) - mp.cosh(t)), 0, 2)
    def integrand(tau):
        if tau < mp.mpf(10)**(-12):
            return 1/(be*mp.sqrt(be*mp.sinh(be)/2))
        g, t = G_stable(tau)
        if t < be*mp.mpf(10)**(-5):
            return G2/2*be*mp.sin(tau)
        return (g - G0)/t**2*be*mp.sin(tau)
    I = mp.quad(integrand, [0, mp.pi/4, mp.pi/2])
    return mp.sqrt(2)/mp.pi*(G0/be - I)
def Phi(m, nu, be):
    return mp.re(mp.cosh(be/2)**(2*m)*mp.hyp2f1(0.5+m+1j*nu, 0.5+m-1j*nu, 1, -mp.sinh(be/2)**2))
def I_direct(n, m, be, kap):
    """regulated moment \int_0^infty e^{-kap nu} nu^n Phi_{m,m} dnu by direct integration of the 2F1 form"""
    numax = 22/kap
    return mp.quad(lambda nu: mp.exp(-kap*nu)*nu**n*Phi(m, nu, be), mp.linspace(0, numax, int(numax*max(1, be)/3) + 2))
def I_MD(n, m, be, kap):
    """the same regulated moment from the Chebyshev-Mehler representation:
    (sqrt2/pi) \int_0^beta G_m(t) Re[n!/(kap - i t)^{n+1}] dt"""
    def integrand(tau):
        if tau < mp.mpf(10)**(-12):
            return mp.re(mp.factorial(n)/(kap - 1j*be)**(n+1))*be/mp.sqrt(be*mp.sinh(be)/2)
        t = be*mp.cos(tau); dt = 2*be*mp.sin(tau/2)**2
        g = T2m(m, mp.cosh(t/2)/mp.cosh(be/2))/mp.sqrt(2*mp.sinh((be + t)/2)*mp.sinh(dt/2))
        return g*mp.re(mp.factorial(n)/(kap - 1j*t)**(n+1))*be*mp.sin(tau)
    brk = [mp.mpf(0)] + [mp.acos(min(mp.mpf(1), c*kap/be)) for c in (30, 10, 3, 1, mp.mpf('0.3'))] + [mp.pi/2]
    brk = sorted(set(brk))
    return mp.sqrt(2)/mp.pi*mp.quad(integrand, brk)
print("\n(b1) regulated moments at kappa = 0.4: direct nu-integral of the 2F1 form vs Chebyshev-Mehler t-integral")
mp.mp.dps = 15
for m in (2, 4):
    for be in (mp.mpf('0.6'), mp.mpf('1.5')):
        for n in (0, 1, 2):
            a = I_direct(n, m, be, mp.mpf('0.4')); b_ = I_MD(n, m, be, mp.mpf('0.4'))
            print(f"   m={m} beta={mp.nstr(be,2)} n={n}:  {mp.nstr(a, 11):>16}  {mp.nstr(b_, 11):>16}")
mp.mp.dps = 20
print("\n(b2) kappa -> 0 (Richardson from kappa = 0.004, 0.002, 0.001 of the t-integral) vs closed forms of Theorem K")
def Q_lim(n, m, be):
    a, b_, c = (I_MD(n, m, be, mp.mpf(k)) for k in ('0.004', '0.002', '0.001'))
    r1, r2 = 2*b_ - a, 2*c - b_
    return (4*r2 - r1)/3
for m in (0, 2, 4):
    for be in map(mp.mpf, ('0.6', '1.0', '2.0')):
        cf = (Q0(m, be), Q1(m, be), Q2(m, be)); dr = tuple(Q_lim(n, m, be) for n in (0, 1, 2))
        print(f"   m={m} beta={mp.nstr(be,2)}:  Q0 {mp.nstr(cf[0],8):>12} [{mp.nstr(dr[0],8):>12}]   Q1 {mp.nstr(cf[1],8):>12} [{mp.nstr(dr[1],8):>12}]   Q2 {mp.nstr(cf[2],8):>12} [{mp.nstr(dr[2],8):>12}]")
print("\n   small-beta limit of Q1^{(m)}: values at beta = 0.02 versus 1/24 - m^2/2")
for m in (0, 2, 4):
    print(f"   m={m}: {mp.nstr(Q1(m, mp.mpf('0.02')), 8)}   {mp.nstr(mp.mpf(1)/24 - mp.mpf(m*m)/2, 8)}")
print("\n(c) table: scalar kernel (pair weights (0,0): both the '++' and '+-' heavy states of the fixed-u MHV sum rule B^(1))")
print("    and Phi_{4,4} (pair weights (4,4): the all-plus sum rule B^(4)); ratios needed for |offdiag| <= diag")
print("    beta      Q1^(0)      Q2^(0)      Q1^(4)      Q2^(4)   |Q1^(4)|/Q1^(0)  |Q2^(4)/Q2^(0)|")
for be in map(mp.mpf, ('0.05', '0.1', '0.2', '0.3', '0.5', '0.8', '1.0', '1.5', '2.0', '3.0')):
    a1, a2, c1, c2 = Q1(0, be), Q2(0, be), Q1(4, be), Q2(4, be)
    print(f"   {mp.nstr(be,3):>5} {mp.nstr(a1,6):>11} {mp.nstr(a2,6):>11} {mp.nstr(c1,6):>11} {mp.nstr(c2,6):>11} {mp.nstr(abs(c1)/a1,6):>12} {mp.nstr(abs(c2/a2),6):>12}")
# weight psi(p) = (1-p)^3 p (65 + 155 p + 47 p^2) of F_g3, arXiv:2201.06602 (3.3a)
import sympy as sp
pp = sp.symbols('p'); psi = sp.expand((1-pp)**3*pp*(65 + 155*pp + 47*pp**2)); print("\n   B_2^(1) weight of F_g3:", psi)
psi4 = sp.expand((1-pp)**3*pp*(99 + 282*pp + 49*pp**2)); print("   B_2^(1) weight of F_g4:", psi4)
a1, a2 = 65, -40
print("   leading curved-regime action a1 Q1/M + a2 Q2/M^2 for (a1,a2) = (65,-40):")
for M in (10, 100, 1000):
    row = []
    for be in map(mp.mpf, ('0.05', '0.2', '0.5', '1.0', '2.0', '3.0')):
        row.append((mp.nstr(a1*Q1(0, be)/M + a2*Q2(0, be)/M**2, 4), mp.nstr(a1*Q1(4, be)/M + a2*Q2(4, be)/M**2, 4)))
    print(f"   M={M}: (scalar kernel, Phi_44 kernel) at beta = 0.05,0.2,0.5,1,2,3:", row)
