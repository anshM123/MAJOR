"""
lsl/lsl_lrinf.py (LARGESPIN_L.md §3, Lemma LR-inf): the leading-Regge per-state action at every impact parameter beta >= beta_0, all
masses m >= Lam, with an explicit remainder.

For k = 2, 3 (psi = psi_k, a1 = psi'(0), a2 = psi''(0), g = psi/p, X = Lam beta, D = cosh(beta) - 1, K0 = D^{-1/2}):
   Omega_k(beta) = a1 M1(beta)/Lam^2 - a2/(2 bhat^3) + E_k(beta),       bhat = 2 Lam sinh(beta/2),
   |E_k| <= (sqrt2/pi) [ T_a + T_b + T_c ],
     T_a = (K0/Lam) d1(X),          d1(X) = sum_{j>=2 even} |g^(j)(0)|/X^{j+1} + sum_{j>=5} |g^(j)(1)|/X^{j+1}
     T_b = K0/(4 D Lam^3) * sum_c,   sum_c = int_X^inf x^2 R_inf(x) dx  (closed form)
     T_c = int_0^{beta/2} R(Lam t) rho_b(t) dt + R_inf(Lam beta/2) (pi/sqrt2) P_{-1/2}(cosh beta)
   with R_inf(x) = sum_{j>=3 odd} |psi^(j)(0)|/x^{j+1} + sum_{j>=5} |psi^(j)(1)|/x^{j+1} >= |r(x)|, r := Phi + a1/x^2 (exact, x > 0),
        R(x) = min(R_inf(x), a1/x^2 + Phi_max),  Phi_max = int|psi|,
        rho_b(t) = K0 [ t^4 cosh(t)/(48 D) + (3/8) 2^{5/2} ((cosh t - 1)/D)^2 ]  >= K(t) - K0 - K0 t^2/(4D)  on t <= beta/2.
Proof: LARGESPIN_L.md §3 (Mehler-Dirichlet; f.p. identity for M1; exact closed form of Phi; Taylor bounds).  The script
  (i) checks the two identities used: (sqrt2/pi)[K0/beta - int_0^beta (K - K0)/t^2] = M1(beta);  C_inf := int_0^inf (x^2 Phi + a1) dx = -(pi/2) a2,
  (ii) evaluates the bound on a beta grid and the ratio |E|/main, (iii) compares with the direct Omega where the direct quadrature converges.
Output: cache/lrinf.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
from scipy.integrate import quad
import lsl_core as C

DG = 1000
L = C.lam(DG)
OUT = os.path.join(C.CACHE, 'lrinf.txt')
lines = []


def pr(s=''):
    print(s, flush=True); lines.append(s)


def coeffs(k):
    a0 = [abs(v) for v in C.d0(k)]; a1 = [abs(v) for v in C.d1(k)]
    # g = psi/p : g^(j)(0) = psi^(j+1)(0)/(j+1);  g^(j)(1) from exact polynomial division
    from fractions import Fraction as Fr
    c = C.psi_coeffs(k)
    gc = c[1:]                                  # psi(0) = 0
    D = len(gc) - 1
    g0 = []; g1 = []
    for j in range(D + 1):
        s0 = Fr(0); s1 = Fr(0)
        for n in range(j, D + 1):
            ff = 1
            for i in range(j):
                ff *= (n - i)
            s1 += gc[n]*ff
            if n == j:
                s0 += gc[n]*ff
        g0.append(abs(float(s0))); g1.append(abs(float(s1)))
    return a0, a1, g0, g1


def R_inf(x, a0, a1):
    x = np.asarray(x, float)
    s = np.zeros_like(x)
    for j in range(3, len(a0), 2):
        s += a0[j]/x**(j + 1)
    for j in range(5, len(a1)):
        s += a1[j]/x**(j + 1)
    return s


def sum_c(X, a0, a1):
    s = 0.0
    for j in range(3, len(a0), 2):
        s += a0[j]*X**(2 - j)/(j - 2)
    for j in range(5, len(a1)):
        s += a1[j]*X**(2 - j)/(j - 2)
    return s


def d1fun(X, g0, g1):
    s = 0.0
    for j in range(2, len(g0), 2):
        s += g0[j]/X**(j + 1)
    for j in range(5, len(g1)):
        s += g1[j]/X**(j + 1)
    return s


def P0(beta):
    return float(mp.legenp(-0.5, 0, mp.cosh(beta), type=3))


def bound(k, beta):
    a0, a1, g0, g1 = coeffs(k)
    A1 = C.d0(k)[1]
    Phimax = float(sum(abs(float(v))/(n + 1) for n, v in enumerate(C.psi_coeffs(k))))   # >= int |psi| (crude, rigorous)
    D = np.cosh(beta) - 1; K0 = D**-0.5; X = L*beta
    Ta = K0/L*d1fun(X, g0, g1)
    Tb = K0/(4*D*L**3)*sum_c(X, a0, a1)
    def integrand(t):
        x = L*t
        R = min(float(R_inf(x, a0, a1)), A1/x**2 + Phimax) if x > 0 else 0.0
        rb = K0*(t**4*np.cosh(t)/(48*D) + 0.375*2**2.5*((np.cosh(t) - 1)/D)**2)
        return R*rb
    pts = sorted(set([0.0, min(beta/2, 5/L), min(beta/2, 30/L), min(beta/2, 300/L), beta/2]))
    Tc1 = sum(quad(integrand, pts[i], pts[i + 1], limit=400, epsabs=0, epsrel=1e-10)[0] for i in range(len(pts) - 1))
    Tc2 = float(R_inf(L*beta/2, a0, a1))*(np.pi/np.sqrt(2))*P0(beta)
    return (np.sqrt(2)/np.pi)*(Ta + Tb + Tc1 + Tc2), dict(Ta=Ta, Tb=Tb, Tc1=Tc1, Tc2=Tc2)


def M1_fp(beta):
    """(sqrt2/pi)[K0/beta - int_0^beta (K - K0)/t^2 dt] via the substitution t = beta(1 - v^2) near the endpoint."""
    mp.mp.dps = 30
    b = mp.mpf(beta); D = mp.cosh(b) - 1; K0 = 1/mp.sqrt(D)
    def g(v):
        if v == 0:
            return mp.mpf(0)
        t = b*(1 - v*v)
        eps = 2*mp.sinh(t/2)**2/D                                   # (cosh t - 1)/D, no cancellation
        if t == 0:
            return K0/(4*D)*2*b*v
        # (1 - eps) = 2 sinh((b+t)/2) sinh((b-t)/2)/D exactly; K - K0 = K0 [ (1-eps)^(-1/2) - 1 ] = K0 expm1(-log1p(-eps)/2)
        one_m = 2*mp.sinh((b + t)/2)*mp.sinh(b*v*v/2)/D
        dK = K0*mp.expm1(-mp.log(one_m)/2) if eps > mp.mpf('0.25') else K0*mp.expm1(-mp.log1p(-eps)/2)
        return dK/t**2*2*b*v
    I = mp.quad(g, [0, 0.1, 0.5, 0.9, 1])
    return float(mp.sqrt(2)/mp.pi*(K0/b - I))


def C_inf_check(k, X):
    """C(X) = int_0^X (x^2 Phi(x) + a1) dx = int_0^1 psi(p) [X^2 sin(pX)/p + 2X cos(pX)/p^2 - 2 sin(pX)/p^3] dp + a1 X."""
    mp.mp.dps = 40
    c = [mp.mpf(v.numerator)/v.denominator for v in C.psi_coeffs(k)]
    psi = lambda p: sum(cn*p**n for n, cn in enumerate(c))
    Xm = mp.mpf(X)
    def f(p):
        if p == 0:
            return mp.mpf(0)
        return psi(p)*(Xm**2*mp.sin(p*Xm)/p + 2*Xm*mp.cos(p*Xm)/p**2 - 2*mp.sin(p*Xm)/p**3)
    I = mp.quad(f, mp.linspace(0, 1, int(X/3) + 3))
    a1 = mp.mpf(C.psi_coeffs(k)[1].numerator)/C.psi_coeffs(k)[1].denominator
    return float(I + a1*Xm)


def main():
    pr(__doc__.split('Output')[0])
    pr("(i) identities")
    for b in [0.1, 1.0, 4.0, 10.0]:
        pr(f"   beta = {b}: M1 (e^(2 pi nu) form) = {C.M1(b):.12e}, f.p. form = {M1_fp(b):.12e}")
    for k in (2, 3):
        a2 = C.d0(k)[2]
        vals = [C_inf_check(k, X) for X in (200.0, 400.0, 800.0)]
        pr(f"   k = {k}: C(X) at X = 200/400/800: {vals[0]:.8f} / {vals[1]:.8f} / {vals[2]:.8f};  -(pi/2) psi''(0) = {-np.pi/2*a2:.8f}")
    pr("\n(ii) bound |E_k|/(a1 M1/Lam^2) on a beta grid (main term a1 M1/Lam^2 - a2/(2 bhat^3); a2 < 0 for k = 2)")
    grid = [0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0, 48.0, 64.0]
    worst = {2: 0, 3: 0}
    for b in grid:
        M1 = C.M1(b) if b <= 16 else None
        if M1 is None:   # h13_cert_m1 region (c) lower bound, beta >= 12
            D = np.cosh(b) - 1; K0 = D**-0.5
            M1 = (np.sqrt(2)/np.pi)*K0*(1/b - 1/(2*np.sinh(b/2)) - np.sqrt(2)*(np.exp(-b/2)/4.7 + 4/b**2) - 2*np.sqrt(2*np.log(2))/(b - np.log(2))**2)
            tag = 'M1 lower bound (h13 region c)'
        else:
            tag = 'M1 numeric'
        bh = 2*L*np.sinh(b/2)
        row = f"   beta = {b:6.3f} (bhat = {bh:9.1f}):"
        for k in (2, 3):
            a1 = C.d0(k)[1]; a2 = C.d0(k)[2]
            E, parts = bound(k, b)
            main = a1*M1/L**2
            rel = E/main
            worst[k] = max(worst[k], rel) if b >= 0.05 else worst[k]
            row += f"  k={k}: |E|/main <= {rel:.2e} [Ta {parts['Ta']/main*np.sqrt(2)/np.pi:.1e} Tb {parts['Tb']/main*np.sqrt(2)/np.pi:.1e} Tc {(parts['Tc1']+parts['Tc2'])/main*np.sqrt(2)/np.pi:.1e}]"
        E2, _ = bound(2, b); E3, _ = bound(3, b)
        A1min = 2*C.d0(2)[1] - C.d0(3)[1]                      # u <= 1, worse sigma
        hyb = (2*E2 + E3)/(A1min*M1/L**2)
        row += f"  | hybrid (u <= 1, worse sigma): (2|E2|+u|E3|)/(A1min M1/Lam^2) <= {hyb:.2e}"
        pr(row + "  " + tag)
    pr(f"   worst ratio on beta >= 0.05: k=2 {worst[2]:.2e}, k=3 {worst[3]:.2e}")
    pr("\n(iii) direct Omega (Mehler-phi quadrature, converged where |conv| small) vs two-term form, relative difference")
    for b in [0.05, 0.1, 0.3, 1.0, 2.0, 4.0, 8.0]:
        row = f"   beta = {b}:"
        for k in (2, 3):
            a1 = C.d0(k)[1]; a2 = C.d0(k)[2]; bh = 2*L*np.sinh(b/2)
            om = C.Omega(k, b, nphi=12000); asym = a1*C.M1(b)/L**2 - a2/(2*bh**3)
            E, _ = bound(k, b)
            row += f"  k={k}: (Omega - form)/form = {(om - asym)/asym:+.2e}, bound/form = {E/asym:.2e}"
        pr(row)
    with open(OUT, 'w') as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == '__main__':
    main()
