"""t7: reality and nu-evenness of C~_3 = Psi~_3/N  (OMEGA_PRIME.md Sec. 8, referee M3) and the d_0 normalisation (M2).

Leading-Regge frame kernel of Psi~_nu: (1,3) antipodal weights w_q(nu) times harmonic F^nu_{q1+q3}; (2,4) kernel A = D5 D1 (modes (+-2,-+2),
mu_R = 0, nu-independent real kernel).  F^nu_mu = P^{|mu|}_{-1/2+i nu}/(1/2+i nu)_{|mu|} and P^mu_{-1/2+i nu} is even in nu, so
F^{-nu}_mu = F^nu_mu (delta)_{|mu|}/(1-delta)_{|mu|}, delta = 1/2 + i nu.  Checks (exact, sympy):
 (i)  shadow:  w_q(-nu) (delta)_{|mu|}/(1-delta)_{|mu|} = [delta/(1-delta)] w_q(nu)   for all odd-mu modes, l = 1, 2
      => Psi~_{-nu} = [delta/(1-delta)] Psi~_nu  at leading Regge order;
 (ii) reality: w_q(-nu) = - conj(w_{-q}(nu)) (real nu) => on real correlators Psi~_{-nu}[G] = -conj(Psi~_nu[G]);
 (iii) N(-nu)/N(nu) = delta/(1-delta) and N(-nu) = -conj N(nu)  (N from t6)
 => C~_{-nu} = C~_nu = conj C~_nu: C~ is real and even in nu (at the level of the leading-Regge kernel).
 (iv) d_0(delta) = (1/pi) int_0^pi sin^{-delta} e^{-i psi} = e^{-i pi/2} 2^delta Gamma(1-delta)/(Gamma((1-delta)/2) Gamma((3-delta)/2)): finite, nonzero on real nu.
"""
import sympy as sp, mpmath as mp
nu = sp.Symbol('nu', real=True)
d = sp.Rational(1, 2) + sp.I*nu
def w(l, q, dd):
    q1, q3 = q; mm = (q1 + q3 - 1)//2; sgn = 1
    if mm < 0: mm = -mm - 1; sgn = -1
    return sp.Integer(-1)**q3*sgn*sp.rf((dd + 1)/2, mm)/sp.rf((3 - dd)/2, mm)
dm = 1 - d                                    # delta at -nu
ok_shadow = ok_conj = True
for l in (1, 2):
    for q1 in range(-l, l+1):
        for q3 in range(-l, l+1):
            if (q1 + q3) % 2 != 1: continue
            mu = abs(q1 + q3)
            lhs = w(l, (q1, q3), dm)*sp.rf(d, mu)/sp.rf(1 - d, mu)
            rhs = d/(1 - d)*w(l, (q1, q3), d)
            ok_shadow &= sp.simplify(sp.expand_func(lhs - rhs)) == 0
            c = sp.simplify(w(l, (q1, q3), dm) + sp.conjugate(w(l, (-q1, -q3), d)))
            ok_conj &= sp.simplify(sp.expand_complex(c)) == 0
print("t7 (i)  shadow identity of the antipodal weights, l = 1, 2, all odd modes:", ok_shadow)
print("t7 (ii) w_q(-nu) = -conj w_{-q}(nu):", ok_conj)
N = lambda x: (2*x + sp.I)*(4*x**2 + 9)/(4*x**2 + 201)
print("t7 (iii) N(-nu)/N(nu) - delta/(1-delta) =", sp.simplify(N(-nu)/N(nu) - d/(1 - d)),
      ";  N(-nu) + conj N(nu) =", sp.simplify(sp.expand_complex(N(-nu) + sp.conjugate(N(nu)))))
print("         l = 1: N1(-nu)/N1(nu) - delta/(1-delta) =", sp.simplify((-2*nu + sp.I)/(2*nu + sp.I) - d/(1 - d)))
mp.mp.dps = 20
print("t7 (iv) d_0: numeric integral vs closed form, and |d_0| on a nu grid:")
for x in [0, 0.5, 1, 2, 4, 8]:
    dl = mp.mpf(1)/2 + 1j*mp.mpf(x)
    num = mp.quad(lambda p: mp.sin(p)**(-dl)*mp.exp(-1j*p), [0, mp.pi/2, mp.pi])/mp.pi
    cf = mp.exp(-1j*mp.pi/2)*2**dl*mp.gamma(1 - dl)/(mp.gamma((1 - dl)/2)*mp.gamma((3 - dl)/2))
    print(f"   nu = {x}: quad = {mp.nstr(num, 10)}, closed = {mp.nstr(cf, 10)}, |d0| = {mp.nstr(abs(cf), 6)}")
print("   zeros of d_0 would need poles of Gamma((1-delta)/2) or Gamma((3-delta)/2): (1-delta)/2 = 1/4 - i nu/2, never a non-positive integer for real nu.")
