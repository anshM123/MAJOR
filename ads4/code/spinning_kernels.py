r"""Regge-channel harmonic kernels for external spinning operators in CFT3 (transverse space H^2, group SL(2,R)).
Phi^nu_{m,n}(beta) := < e^{i m theta}, pi_nu(a_beta) e^{i n theta} >, pi_nu the spherical principal series on L^2(S^1):
   Phi^nu_{m,n}(beta) = (1/2pi) \int_0^{2pi} e^{-i m theta} (cosh b - sinh b cos theta)^{-1/2 - i nu} e^{i n theta'} d theta,
   e^{i theta'} = (e^{i theta} cosh(b/2) - sinh(b/2)) / (cosh(b/2) - e^{i theta} sinh(b/2)).
Checks: (1) Phi_{0,0} = P_{-1/2+i nu}(cosh b);  (2) Phi_{m,m} = cosh(b/2)^{2m} 2F1(1/2+m+i nu, 1/2+m-i nu; 1; -sinh^2(b/2));
(3) flat limit Phi_{m,m} -> J_0(nu b), Phi_{0,m} ~ J_m;  (4) Mehler-Dirichlet form
   Phi_{m,m} = (sqrt2/pi) \int_0^b cos(nu t) c_m(t) (cosh b - cosh t)^{-1/2} dt,  c_m(t) = cos(m chi(t));
(5) zeroth moment  Q_0^{(m,m)}(b) = \int_0^infty Phi^nu_{m,m} d nu = cos(2 m gd(b/2)) / (2 sinh(b/2)),  gd = Gudermannian."""
import mpmath as mp
mp.mp.dps = 20
def Phi_int(m, n, nu, b):
    c, s = mp.cosh(b/2), mp.sinh(b/2)
    def f(th):
        e = mp.exp(1j*th); ep = (e*c - s)/(c - e*s)
        return mp.exp(-1j*m*th)*(mp.cosh(b)-mp.sinh(b)*mp.cos(th))**(-0.5-1j*nu)*ep**n
    return mp.quad(f, mp.linspace(0, 2*mp.pi, 9))/(2*mp.pi)
def Phi_diag(m, nu, b):
    return mp.re(mp.cosh(b/2)**(2*m)*mp.hyp2f1(0.5+m+1j*nu, 0.5+m-1j*nu, 1, -mp.sinh(b/2)**2))
def chi(t, b):
    c, s = mp.cosh(b/2), mp.sinh(b/2); ct = (mp.cosh(b)-mp.exp(t))/mp.sinh(b); st = mp.sqrt(max(0, 1-ct*ct))
    return 2*mp.atan2(s*st, c - s*ct)
def Phi_md(m, nu, b):
    return (mp.sqrt(2)/mp.pi)*mp.quad(lambda t: mp.cos(nu*t)*mp.cos(m*chi(t, b))/mp.sqrt(mp.cosh(b)-mp.cosh(t)), [0, b/2, b])
print("(1),(2),(4): integral definition vs closed form vs Mehler-Dirichlet form")
for (m, nu, b) in ((0, 1.3, 0.8), (4, 1.3, 0.8), (4, 3.0, 2.0), (2, 0.7, 1.5)):
    nu, b = mp.mpf(nu), mp.mpf(b); a = Phi_int(m, m, nu, b)
    ref = mp.re(mp.legenp(-0.5+1j*nu, 0, mp.cosh(b), type=3)) if m == 0 else Phi_diag(m, nu, b)
    print(f"  m={m} nu={mp.nstr(nu,3)} b={mp.nstr(b,3)}:  integral {mp.nstr(a, 12)}   closed form {mp.nstr(ref, 12)}   Mehler-Dirichlet {mp.nstr(Phi_md(m, nu, b), 12)}")
print("chi even in t?", [mp.nstr(chi(mp.mpf(t), mp.mpf(1.2)) - chi(-mp.mpf(t), mp.mpf(1.2)), 5) for t in (0.3, 0.9)])
print("(3) flat limit, m=4: nu=200, b=0.02 (nu b = 4):", mp.nstr(Phi_diag(4, mp.mpf(200), mp.mpf('0.02')), 8), " J_0(4) =", mp.nstr(mp.besselj(0, 4), 8),
      "  P_{-1/2+i nu} =", mp.nstr(mp.re(mp.legenp(-0.5+200j, 0, mp.cosh(mp.mpf('0.02')), type=3)), 8))
print("(5) zeroth moment of the diagonal kernels: Abel-regulated numerical integral vs cos(2 m gd(b/2))/(2 sinh(b/2))")
gd = lambda x: 2*mp.atan(mp.tanh(x/2))
for m in (0, 2, 4):
    for b in map(mp.mpf, ('0.3', '0.6', '1.0', '2.0')):
        vals = []
        for kap in (mp.mpf('0.08'), mp.mpf('0.04')):
            vals.append(mp.quad(lambda nu: mp.exp(-kap*nu)*Phi_diag(m, nu, b), mp.linspace(0, 240, 241)))
        num = 2*vals[1]-vals[0]                                  # Richardson in kappa
        exact = mp.cos(2*m*gd(b/2))/(2*mp.sinh(b/2))
        print(f"  m={m} beta={mp.nstr(b,3):>4}:  numerical {mp.nstr(num, 7):>12}   formula {mp.nstr(exact, 7):>12}")
print("sign changes of Q_0^{(4,4)}: zeros at gd(beta/2) = (2j+1) pi/16, i.e. beta =", [mp.nstr(2*mp.asinh(mp.tan((2*j+1)*mp.pi/16)), 5) for j in range(4)])
