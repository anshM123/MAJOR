"""Step 7: extreme-component norm N_{-J} for symbolic J, and assembly of the absolute kernels.
(a) N_{-J} = pi^{5/2} 2^{4-2Delta-J} (Delta-1)_J / (Gamma(Delta+J) Gamma(Delta-1/2)): the Fourier transform of
    (zb.zk x^2 - 2 zb.x zk.x)^J (x^2)^{-Delta-J}, zb = e_+/2, zk = e_-/2, at p = (1,0,0) reduces (zeta.p = 0, null zeta) to
    pi^{5/2} 2^{4-2Delta-J}/Gamma(Delta-1/2) sum_k C(J,k)(-1)^k k!/Gamma(Delta+k) = ... 2F1(-J,1;Delta;1)/Gamma(Delta) (Gauss).
(b) |K_0|^2 = 12 |M_{-J}|^2 / N_{-J} / W0_J  and  |K_4|^2 = 12 |M_J|^2 / N_J / W4_{J-4}  reproduce the KERNELS_D3 closed
    forms (symbolic check with gammasimp at symbolic Delta and J)."""
from common import *
D = sp.Symbol('Delta', positive=True)
Jp = sp.Symbol('J', positive=True, integer=True)
k = sp.Symbol('k', integer=True, nonnegative=True)

# (a) the sum, symbolically (sympy evaluates the terminating 2F1 by Gauss/Chu-Vandermonde)
ssum = sp.summation(sp.binomial(Jp, k)*(-1)**k*sp.factorial(k)/sp.gamma(D + k), (k, 0, Jp))
target = sp.rf(D - 1, Jp)/sp.gamma(D + Jp)
print("(a) sum_k C(J,k)(-1)^k k!/Gamma(Delta+k) - (Delta-1)_J/Gamma(Delta+J) =",
      sp.simplify(sp.gammasimp(sp.hyperexpand(ssum) - target)))
for Jv in (4, 7, 20):
    print(f"    J={Jv}: exact check", sp.simplify(sum(sp.binomial(Jv, q)*(-1)**q*sp.factorial(q)/sp.gamma(D + q) for q in range(Jv + 1)) - target.subs(Jp, Jv)) == 0)

# (b) assembly
s_ = (D + Jp)/2; h_ = (D - Jp)/2
pi = sp.pi
def N(m):   # KERNELS_D3 3.3 (used for the ladder ratio); at m = -J it is (a)
    return pi**sp.Rational(5, 2)*2**(4 - 2*D - Jp)*sp.binomial(2*Jp, Jp + m)*sp.rf(D - Jp - 1, Jp - abs(m))*sp.rf(D - 1, abs(m))/(sp.gamma(D + Jp)*sp.gamma(D - sp.Rational(1, 2)))
NmJ = pi**sp.Rational(5, 2)*2**(4 - 2*D - Jp)*sp.rf(D - 1, Jp)/(sp.gamma(D + Jp)*sp.gamma(D - sp.Rational(1, 2)))
W0J = sp.binomial(Jp + 3, 3)**2*sp.rf(1, 2*Jp)/sp.rf(D - 1 - Jp, 2*Jp)            # W0_k at k = J
W4 = sp.binomial(Jp - 1, 3)**2*sp.rf(1, Jp - 4)*sp.rf(Jp + 5, Jp - 4)/(sp.rf(D - 1 - Jp, Jp - 4)*sp.rf(D + 3, Jp - 4))  # W4_k at k = J-4
# |M_{-J}|, |M_J| from ap5 (moduli; e^{i pi s}, (-1)^{J/2} have modulus 1 for real s)
absM = {'pp': 4*pi**3*2**(-Jp)*sp.factorial(Jp + 3)/(sp.gamma(s_)*sp.gamma(D + 3)*sp.gamma(1 - h_)),
        'odd': 8*pi**3*2**(-Jp)*sp.factorial(Jp + 3)/(sp.gamma(s_ + sp.Rational(1, 2))*sp.gamma(D + 3)*sp.gamma(sp.Rational(3, 2) - h_)),
        'pm': 2*pi**3*2**(-Jp)*sp.factorial(Jp - 1)/(sp.gamma(s_)*sp.gamma(D - 1)*sp.gamma(1 - h_)),
        'oddJ': 8*pi**3*2**(-Jp)*sp.factorial(Jp - 1)/(sp.gamma(s_ + sp.Rational(1, 2))*sp.gamma(D - 1)*sp.gamma(sp.Rational(3, 2) - h_))}
closed = {   # KERNELS_D3 5.3
    'pp': 27*pi**sp.Rational(7, 2)*2**(4 + 2*D - Jp)*sp.rf(D - Jp - 1, Jp)*sp.gamma(D + Jp)*sp.gamma(D - sp.Rational(1, 2))/(sp.binomial(2*Jp, Jp)*sp.gamma(s_)**2*sp.gamma(D + 3)**2*sp.gamma(1 - h_)**2),
    'odd': 27*pi**sp.Rational(7, 2)*2**(6 + 2*D - Jp)*sp.rf(D - Jp - 1, Jp)*sp.gamma(D + Jp)*sp.gamma(D - sp.Rational(1, 2))/(sp.binomial(2*Jp, Jp)*sp.gamma(s_ + sp.Rational(1, 2))**2*sp.gamma(D + 3)**2*sp.gamma(sp.Rational(3, 2) - h_)**2),
    'pm': 27*pi**sp.Rational(7, 2)*2**(2 + 2*D - Jp)*sp.rf(D - Jp - 1, Jp - 4)*sp.gamma(D + Jp)*sp.gamma(D - sp.Rational(1, 2))/(sp.binomial(2*Jp, Jp + 4)*sp.rf(D - 1, 4)*sp.gamma(s_)**2*sp.gamma(D - 1)**2*sp.gamma(1 - h_)**2),
    'oddJ': 27*pi**sp.Rational(7, 2)*2**(6 + 2*D - Jp)*sp.rf(D - Jp - 1, Jp - 4)*sp.gamma(D + Jp)*sp.gamma(D - sp.Rational(1, 2))/(sp.binomial(2*Jp, Jp + 4)*sp.rf(D - 1, 4)*sp.gamma(s_ + sp.Rational(1, 2))**2*sp.gamma(D - 1)**2*sp.gamma(sp.Rational(3, 2) - h_)**2)}
print("(a') N_{-J} formula = KERNELS_D3 N_{O,m} at m = -J:", sp.simplify(sp.gammasimp(N(-Jp)/NmJ)) == 1)
for name in ('pp', 'odd', 'pm', 'oddJ'):
    if name in ('pp', 'odd'):
        K2 = 12*absM[name]**2/NmJ/W0J
    else:
        K2 = 12*absM[name]**2/NmJ/W4          # N_{+J} = N_{-J}
    r = sp.gammasimp(sp.combsimp(sp.expand_func((K2/closed[name]).rewrite(sp.gamma))))
    print(f"(b) {name:5s}: [12|M_ext|^2/N_ext/W_ext] / [KERNELS_D3 closed form] = {sp.simplify(r)}")
