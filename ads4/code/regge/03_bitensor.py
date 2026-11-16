r"""Check 3: both pairs spinning.  Bi-tensor derivatives of the scalar harmonic function on H^2 and the
SL(2,R) principal-series matrix coefficients.

Poincare disk, ds^2 = e^{2 sigma} dw dwbar, e^{2 sigma} = 4/(1-|w|^2)^2.  For a function f the helicity-m
derivative is the (w...w) component of nabla^m f:  nabla_w T_{w^k} = (d_w - 2k d_w sigma) T_{w^k},
2 d_w sigma = 2 wbar/(1-|w|^2), and the unit null vectors are eps_- = sqrt2 e^{-sigma} d_w, eps_+ = sqrt2 e^{-sigma} d_wbar.

(A) Raising lemma (symbolic, m <= 5): for the plane wave E_s(w; b) = ((1-|w|^2)/|1 - w bbar|^2)^s, |b| = 1,
        nabla_w^m E_s = (s)_m A^m E_s,   A = (bbar - wbar)/((1-|w|^2)(1 - w bbar)).
(B) Two-point statement (symbolic derivatives, numerical evaluation, m = 1..4): with F = P_{-1/2+i nu}(cosh L(x,y)),
    x = 0, y = tanh(beta/2) on the real axis (frames parallel along the geodesic),
        (eps_-.nabla_x)^m (eps_+.nabla_y)^m F = 2^{-m} prod_{j<m} (nu^2+(j+1/2)^2) Phi^nu_{m, m}(beta),
        (eps_-.nabla_x)^m (eps_-.nabla_y)^m F = 2^{-m} prod_{j<m} (nu^2+(j+1/2)^2) Phi^nu_{m,-m}(beta),
        (eps_-.nabla_x)^m F               = 2^{-m/2} (1/2+i nu)_m Phi^{-nu}_{m,0}(beta),
    where Phi^nu_{m,n}(beta) = (1/2pi) \int e^{-i m th} (cosh b - sinh b cos th)^{-1/2-i nu} e^{i n th'} dth is the
    matrix coefficient of the unitary principal series (definition of ads4/SPINNING_KERNELS.md)."""
import sympy as sp, mpmath as mp
mp.mp.dps = 30
# ---------- (A) raising lemma ----------
w, wb, bb, s = sp.symbols('w wb bb s')          # bb = conj(b), b = 1/bb on the unit circle
E = ((1 - w*wb)/((1 - w*bb)*(1 - wb/bb)))**s
A = (bb - wb)/((1 - w*wb)*(1 - w*bb))
T = E; ok = True
for m in range(1, 6):
    T = sp.diff(T, w) - 2*(m-1)*wb/(1 - w*wb)*T
    ratio = sp.simplify(T/(sp.rf(s, m)*A**m*E))
    ok = ok and (ratio == 1)
    print(f"(A) m={m}: nabla_w^m E_s / ((s)_m A^m E_s) =", ratio)
# ---------- (B) two-point derivatives ----------
w1, w2, w2b = sp.symbols('w1 w2 w2b')
g = sp.symbols('g0:10')                         # g_n = d^n P/du^n at u = cosh L
u = 1 - 2*(w1 - w2)*w2b/(1 - w2*w2b)            # cosh L(w1, w2) at w1bar = 0 (no d/dw1bar is taken)
def D(expr, var, uu):                           # derivative of sum_n c_n(w) g_n(uu)
    out = sp.diff(expr, var)
    for n in range(9):
        out += sp.diff(expr, g[n])*g[n+1]*sp.diff(uu, var)
    return out
u0 = u.subs(w1, 0)
def two_point(m, kind):
    T = g[0]
    for k in range(m):                          # at x: wbar_1 = 0 so nabla_w = d_w exactly
        T = D(T, w1, u)
    T = T.subs(w1, 0)
    for k in range(m):                          # at y
        if kind == 'same':                      # nabla_w at y  -> Phi_{m,-m}
            T = D(T, w2, u0) - 2*k*w2b/(1 - w2*w2b)*T
        elif kind == 'opp':                     # nabla_wbar at y -> Phi_{m,m}
            T = D(T, w2b, u0) - 2*k*w2/(1 - w2*w2b)*T
    return sp.simplify(T)
def gnum(n, nu, x):                             # d^n/dx^n P_{-1/2+i nu}(x)
    a, b = mp.mpf(1)/2 + 1j*nu, mp.mpf(1)/2 - 1j*nu
    return mp.rf(a, n)*mp.rf(b, n)/mp.factorial(n)*(-mp.mpf(1)/2)**n*mp.hyp2f1(a+n, b+n, 1+n, (1-x)/2)
def Phi_int(m, n, nu, b):                       # unitary matrix coefficient, integral definition
    c, sh = mp.cosh(b/2), mp.sinh(b/2)
    def f(th):
        e = mp.exp(1j*th); ep = (e*c - sh)/(c - e*sh)
        return mp.exp(-1j*m*th)*(mp.cosh(b) - mp.sinh(b)*mp.cos(th))**(-0.5 - 1j*nu)*ep**n
    return mp.quad(f, mp.linspace(0, 2*mp.pi, 17))/(2*mp.pi)
def Phi_diag(m, nu, b):
    return mp.cosh(b/2)**(2*m)*mp.hyp2f1(0.5+m+1j*nu, 0.5+m-1j*nu, 1, -mp.sinh(b/2)**2)
print("\n(B) derivative kernels versus matrix coefficients (columns: derivative, prefactor*Phi, ratio-1)")
worst = 0
for m in (1, 2, 3, 4):
    exprs = {'opp': two_point(m, 'opp'), 'same': two_point(m, 'same'), 'one': two_point(m, 'none')}
    for (nu, beta) in ((mp.mpf('0.7'), mp.mpf('0.6')), (mp.mpf('1.3'), mp.mpf('1.5')), (mp.mpf('3.1'), mp.mpf('2.2'))):
        t = mp.tanh(beta/2); x = mp.cosh(beta)
        sub = {w2: t, w2b: t}; sub.update({g[n]: gnum(n, nu, x) for n in range(10)})
        def ev(e):
            z = sp.N(e.subs({k: sp.Float(mp.nstr(mp.re(v), 30), 30) + sp.I*sp.Float(mp.nstr(mp.im(v), 30), 30) for k, v in sub.items()}), 28)
            return mp.mpc(str(sp.re(z)), str(sp.im(z)))
        fx, fy = (mp.sqrt(2)/2)**m, (mp.sqrt(2)*(1 - t*t)/2)**m     # (sqrt2 e^{-sigma})^m at x = 0 and at y = t
        poly = mp.mpf(1)
        for j in range(m):
            poly *= nu**2 + (j + mp.mpf(1)/2)**2
        d_opp = fx*fy*ev(exprs['opp']); r_opp = poly/2**m*Phi_diag(m, nu, beta)
        d_same = fx*fy*ev(exprs['same']); r_same = poly/2**m*Phi_int(m, -m, nu, beta)
        d_one = fx*ev(exprs['one']); r_one = mp.rf(mp.mpf(1)/2 + 1j*nu, m)/mp.sqrt(2)**m*Phi_int(m, 0, -nu, beta)
        p_one = mp.legenp(-0.5+1j*nu, m, x, type=3)/mp.sqrt(2)**m
        for nm, a, b in (("(m, m)", d_opp, r_opp), ("(m,-m)", d_same, r_same), ("(m, 0)", d_one, r_one)):
            worst = max(worst, abs(a - b)/abs(b))
            print(f"  m={m} nu={mp.nstr(nu,3)} beta={mp.nstr(beta,3)} {nm}: {mp.nstr(mp.re(a),14):>20} {mp.nstr(mp.re(b),14):>20}  |rel.dev.|={mp.nstr(abs(a-b)/abs(b),3)}")
        worst = max(worst, abs(d_one - (-1)**m*p_one)/abs(p_one))
print("max relative deviation over all cases:", mp.nstr(worst, 5), "  [(m,0) also equals (-1)^m 2^{-m/2} P^m_{-1/2+i nu}(cosh beta)]")
# imaginary parts of Phi_{m,-m} and evenness in nu
for m in (2, 4):
    a = Phi_int(m, -m, mp.mpf('1.3'), mp.mpf('1.5')); b = Phi_int(m, -m, -mp.mpf('1.3'), mp.mpf('1.5'))
    print(f"Phi_({m},-{m}): value {mp.nstr(a, 12)}; at -nu {mp.nstr(b, 12)}")
# flat limits
nu, beta = mp.mpf(300), mp.mpf('0.02')
print("flat limit nu*beta = 6:  Phi_44 =", mp.nstr(mp.re(Phi_diag(4, nu, beta)), 8), " J_0(6) =", mp.nstr(mp.besselj(0, 6), 8),
      ";  Phi_4,-4 =", mp.nstr(mp.re(Phi_int(4, -4, nu, beta)), 8), " J_8(6) =", mp.nstr(mp.besselj(8, 6), 8),
      ";  Phi_4,0 =", mp.nstr(Phi_int(4, 0, nu, beta), 8), " J_4(6) =", mp.nstr(mp.besselj(4, 6), 8))
