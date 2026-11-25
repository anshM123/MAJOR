"""[certified-exact check] author's dressing (ANGULAR_ROUTE §45) vs Phi2 = w^3 sum_k V_{b,k} g^{(3)}_k, in exact rationals at rational nu^2
(K_3(nu+i) is a polynomial in (nu+i)^2 -> Gaussian rationals; we check the identity as polynomials via exact evaluation)."""
import sys, os, random
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
def gco(p, s, K):   # s = nu^2 + 1/4 (may be complex rational given as (re, im) Fractions) -> use python complex of Fractions via pairs
    g = [(Fr(1), Fr(0))]; gm = (Fr(0), Fr(0))
    for k in range(K):
        c = (2*k*k + 2*p*k + p - s[0], -s[1])
        prod = (c[0]*g[k][0] - c[1]*g[k][1], c[0]*g[k][1] + c[1]*g[k][0])
        prev = g[k - 1] if k else gm
        q = (k + p - 1)**2
        g.append(((prod[0] - q*prev[0])/(k + 1)**2, (prod[1] - q*prev[1])/(k + 1)**2))
    return g
def Kp(p, s, U):
    g = gco(p, s, len(U))
    return (sum(U[k]*g[k][0] for k in range(len(U))), sum(U[k]*g[k][1] for k in range(len(U))))
random.seed(7); bad = 0; tot = 0
for _ in range(40):
    b = random.choice([0, 4]); J = b + random.randint(1, 14); Delta = J + random.randint(3, 40)
    T = (Delta - 1)*(Delta - 2); C = J*(J + 1); n = J - b
    U = [Fr(1)]
    for k in range(n): u = (b + k)*(b + k + 1); U.append(U[-1]*Fr(C - u, T - u))
    U.append(Fr(0))
    x = Fr(random.randint(0, 500), random.randint(1, 9))      # nu^2
    # Phi2 / w^3 = sum_k V_k g3_k
    V = []
    for k in range(n + 2):
        um = (b - 1)*b if k == 0 else (b + k - 1)*(b + k)
        v = U[k]*(1 + (Fr((k + 3)**2, 2) - um)/T)
        if k >= 1: v += Fr(k*k, 2*T)*U[k - 1]
        V.append(v)
    g3 = gco(3, (x + Fr(1, 4), Fr(0)), n + 2)
    lhs = sum(V[k]*g3[k][0] for k in range(n + 2))
    K3 = Kp(3, (x + Fr(1, 4), Fr(0)), U)[0]; K4 = Kp(4, (x + Fr(1, 4), Fr(0)), U)[0]
    # K3(nu+i): s = (nu+i)^2 + 1/4 = x - 1 + 1/4 + 2 i nu  -> needs nu itself; choose x a perfect square rational
    nu = Fr(random.randint(0, 40), random.randint(1, 5)); x = nu*nu
    g3 = gco(3, (x + Fr(1, 4), Fr(0)), n + 2); lhs = sum(V[k]*g3[k][0] for k in range(n + 2))
    K3 = Kp(3, (x + Fr(1, 4), Fr(0)), U)[0]; K4 = Kp(4, (x + Fr(1, 4), Fr(0)), U)[0]
    K3s = Kp(3, (x - 1 + Fr(1, 4), 2*nu), U)          # K3(nu + i), complex
    l = (Fr(-1, 2), nu)                                 # ell = -1/2 + i nu
    # (7-2b)/T [ (3+l) K4 - (3 + l/2) K3 - (l/2) K3(nu+i) ]
    t1 = ((3 + l[0])*K4, l[1]*K4)
    t2 = ((3 + l[0]/2)*K3, (l[1]/2)*K3)
    t3 = ((l[0]*K3s[0] - l[1]*K3s[1])/2, (l[0]*K3s[1] + l[1]*K3s[0])/2)
    br = (t1[0] - t2[0] - t3[0], t1[1] - t2[1] - t3[1])
    rhs_re = (1 + Fr(9 - b*(b - 1), T))*K3 - (x + Fr(25, 4))/(2*T)*K4 + Fr(7 - 2*b, 1)/T*br[0]
    rhs_im = Fr(7 - 2*b, 1)/T*br[1]
    tot += 1
    if lhs != rhs_re: bad += 1
    if _ < 3: print('b', b, 'J', J, 'Delta', Delta, 'nu', nu, ' lhs-rhs_re =', float(lhs - rhs_re), ' rhs_im =', float(rhs_im), ' lhs=', float(lhs))
print('dressing identity (real part): mismatches', bad, '/', tot)
