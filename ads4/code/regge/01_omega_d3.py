r"""Check 1: the transverse harmonic function of conformal Regge theory in d = 3.
Costa-Hansen-Penedones arXiv:1707.07689, (A.13), (A.17), h = d/2:
    Pi_{Delta-h}(rho) = pi^{1-h} Gamma(Delta-1) / (2 Gamma(Delta-h+1)) e^{(1-Delta) rho} 2F1(h-1, Delta-1; Delta-h+1; e^{-2 rho}),
    Omega_{i nu}(rho) = (i nu / 2 pi) (Pi_{i nu} - Pi_{-i nu}),   Delta = h + i nu.
Kulaxizi-Parnachev-Zhiboedov arXiv:1705.02934 (A.1) (with sinh(pi nu)):
    Omega = nu sinh(pi nu) Gamma(h-1+i nu) Gamma(h-1-i nu) / (2^{d-1} pi^{(d+1)/2} Gamma((d-1)/2)) 2F1(h-1+i nu, h-1-i nu; h-1/2; -sinh^2(rho/2)).
Claim (d = 3):  Omega_{i nu}(L) = nu tanh(pi nu) / (4 pi) * P_{-1/2+i nu}(cosh L),
and (Laplacian_{H^2} + nu^2 + 1/4) Omega = 0."""
import mpmath as mp
mp.mp.dps = 30
h = mp.mpf(3)/2
def Pi(a, rho):          # a = Delta - h
    D = h + a
    return mp.pi**(1-h)*mp.gamma(D-1)/(2*mp.gamma(D-h+1))*mp.exp((1-D)*rho)*mp.hyp2f1(h-1, D-1, D-h+1, mp.exp(-2*rho))
def Omega_CHP(nu, rho):
    return 1j*nu/(2*mp.pi)*(Pi(1j*nu, rho) - Pi(-1j*nu, rho))
def Omega_KPZ(nu, rho, d=3):
    hh = mp.mpf(d)/2
    return nu*mp.sinh(mp.pi*nu)*mp.gamma(hh-1+1j*nu)*mp.gamma(hh-1-1j*nu)/(2**(d-1)*mp.pi**(mp.mpf(d+1)/2)*mp.gamma(mp.mpf(d-1)/2)) \
        * mp.hyp2f1(hh-1+1j*nu, hh-1-1j*nu, hh-mp.mpf(1)/2, -mp.sinh(rho/2)**2)
def P(nu, x):
    return mp.hyp2f1(mp.mpf(1)/2+1j*nu, mp.mpf(1)/2-1j*nu, 1, (1-x)/2)
worst = 0
for nu in map(mp.mpf, ('0.3', '1.7', '4.2')):
    for L in map(mp.mpf, ('0.2', '1.1', '2.5')):
        a = Omega_CHP(nu, L); b = Omega_KPZ(nu, L); c = nu*mp.tanh(mp.pi*nu)/(4*mp.pi)*P(nu, mp.cosh(L))
        worst = max(worst, abs(a-c)/abs(c), abs(b-c)/abs(c))
        print(f"nu={mp.nstr(nu,3)} L={mp.nstr(L,3)}  CHP {mp.nstr(mp.re(a),15)}  KPZ {mp.nstr(mp.re(b),15)}  nu tanh(pi nu)/(4 pi) P {mp.nstr(mp.re(c),15)}")
print("max relative deviation:", mp.nstr(worst, 5))
# eigenvalue equation
nu, L = mp.mpf('1.3'), mp.mpf('0.9')
F = lambda x: mp.re(Omega_CHP(nu, x))
lap = mp.diff(F, L, 2) + mp.coth(L)*mp.diff(F, L)
print("(Lap + nu^2 + 1/4) Omega =", mp.nstr(lap + (nu**2 + mp.mpf(1)/4)*F(L), 5))
