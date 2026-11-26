# r7_k3light.py -- k=3 light normalisation of MIXED_PROFILE_THEOREM: coefficients of (1-z)^{-4}P_{-1/2+i nu}((1+z)/(1-z))
# (fp-cdh at d=1: sum D_k^2 q_k z^k) and the identity 4 g_{2,4} = nu^4-(35/2)nu^2+569/16 (K_B2A4 sec 3 light polynomial),
# hence K_{k=3} = 2 int psi3 [p^4-(35/2)p^2/L^2+(569/16)/L^4] = (2/L^6) int E_3 (4 g_{2,4}) = 8 G_2/L^6 when G_0=G_1=0.
import sympy as sp
z, nu = sp.symbols('z nu')
a = sp.Rational(1, 2)-sp.I*nu
ser = sp.series((1-z)**(a-4)*sum(sp.rf(a, j)**2/sp.factorial(j)**2*z**j for j in range(4)), z, 0, 4).removeO()
g = [sp.expand(sp.simplify(ser.coeff(z, j))) for j in range(4)]
for j in range(4): print('g_%d,4 =' % j, g[j])
# compare with D_k^2 q_k from the paper's recurrence (n+4)^2 q_{n+1} = -(s+4-2(n+2)^2) q_n - n^2 q_{n-1}
s = nu**2+sp.Rational(1, 4); q = [sp.Integer(1)]
for n in range(3):
    q.append(sp.expand((-(s+4-2*(n+2)**2)*q[n]-(n**2*q[n-1] if n >= 1 else 0))/(n+4)**2))
print('g_{k,4} = D_k^2 q_k (fp-cdh d=1):', all(sp.expand(g[k]-sp.binomial(k+3, 3)**2*q[k]) == 0 for k in range(4)))
print('4 g_{2,4} - (nu^4-35/2 nu^2+569/16) =', sp.expand(4*g[2]-(nu**4-sp.Rational(35, 2)*nu**2+sp.Rational(569, 16))))
