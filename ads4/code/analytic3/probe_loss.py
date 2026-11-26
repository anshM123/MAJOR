"""DV: where does K_cell lose against the exact termwise sum?  exact x_j = S^(5/2 - j + i nu) (probe3 method) vs Lemma B with exact
M0, M1 and with the certified M0, M1 (rt2_kernel.PP_box) over a small a-box."""
import sys, math
sys.path.insert(0, '../analytic2'); sys.path.insert(0, '../analytic')
import numpy as np
from common2 import n_of_L
from fastS import S_vec
from flint import arb, acb, ctx
import bchain as BC, rt2_kernel as RK, rt_cert as RC
ctx.prec = 300
L = float(sys.argv[1]); nu = float(sys.argv[2]); a = float(sys.argv[3]); JJ = int(sys.argv[4]) if len(sys.argv) > 4 else 110
b = 0; s = 1
n = n_of_L(L, a, b)
N = arb(n) + arb(a); D = N*(N + s); d = 1 - arb(n*(n + s))/D
I = acb(0, 1); al = acb(2.5, nu)
S0, lg = S_vec(n, a, s, [complex(2.5, nu), complex(3.5, nu)])
Sv = [acb(complex(S0[i]))*arb(lg[i]).exp() for i in range(2)]
Sd = {1: Sv[1], 0: Sv[0]}
for j in range(0, JJ):
    beta = al - j
    Sd[-(j + 1)] = ((d*D + beta*(s - 2*beta - 1))*Sd[-j] + beta*(beta + 1)*Sd[-j + 1])/(D - beta*(beta - s))
xs = [float(abs(Sd[1]*d**(al + 1)).mid())] + [float(abs(Sd[-j]*d**(al - j)).mid()) for j in range(JJ)]
df = float(d.mid())
cs = RC.cj_list(acb(nu), JJ)
B = float(abs((2*I*nu).gamma()*((acb(0.5) + I*nu).rgamma())**2).mid())
c = acb(1); a0 = acb(0.5) - I*nu; b0 = acb(1) - 2*I*nu; tot = acb(0)
for j in range(JJ):
    tot += c*Sd[-j]; c = c*(a0 + j)*(a0 + j)/((b0 + j)*(j + 1))
Btrue = (2*I*nu).gamma()*((acb(0.5) + I*nu).rgamma())**2
print(f'L={L} nu={nu} a={a} n={n} d={df:.5f} kappa={nu/a:.3f}  K_true(J={JJ}) = {float((2*(Btrue*tot).real).mid()):.4g}')
term = [cs[j]*df**(j - 2.5)*xs[j + 1] for j in range(JJ)]
print(' exact termwise |B| sum |c_j| d^{j-5/2}|S^_j| =', B*sum(term), ' peak j', int(np.argmax(term)))
dl, dh = df*(1 - 1e-6), df*(1 + 1e-6)
chB = BC.chain_bounds_multi(2.5, nu, nu + 1e-6, xs[1], xs[0], JJ - 1, dl, dh, a, a + 1e-6)
tB = [cs[j]*df**(j - 2.5)*chB[j + 1] for j in range(JJ - 1)]
print(' Lemma B with exact M0,M1:', B*sum(tB))
M0 = RK.PP_box(2.5, nu, nu + 1.1, a, a + 1, dl, dh); M1 = RK.PP_box(3.5, nu, nu + 1.1, a, a + 1, dl, dh)
print(f' certified M0 = {M0:.3g} (exact {xs[1]:.3g}), M1 = {M1:.3g} (exact {xs[0]:.3g})')
chC = BC.chain_bounds_multi(2.5, nu, nu + 1.1, M0, M1, JJ - 1, dl, dh, a, a + 1)
tC = [cs[j]*df**(j - 2.5)*chC[j + 1] for j in range(JJ - 1)]
print(' Lemma B with certified M0, M1, box:', B*sum(tC))
for j in (0, 2, 5, 10, 20, 30, 40):
    print(f'   j={j}: exact {xs[j+1]:.3g}  B(exactM) {chB[j+1]:.3g}  B(cert) {chC[j+1]:.3g}  |c_j|d^j {cs[j]*df**j:.3g}')
