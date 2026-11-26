"""DV: midpoint solution frames for the beta-chain (|x_j| <= |P_j| M0 + |Q_j| M1 at a point) vs Lemma B multi and exact."""
import sys, math
sys.path.insert(0, '../analytic2'); sys.path.insert(0, '../analytic')
from common2 import n_of_L
from fastS import S_vec
from flint import arb, acb, ctx
import bchain as BC
ctx.prec = 300
L = float(sys.argv[1]); nu = float(sys.argv[2]); a = float(sys.argv[3]); JJ = 41
s = 1; n = n_of_L(L, a, 0)
N = arb(n) + arb(a); D = N*(N + s); d = 1 - arb(n*(n + s))/D
al = acb(2.5, nu)
S0, lg = S_vec(n, a, s, [complex(2.5, nu), complex(3.5, nu)])
Sv = [acb(complex(S0[i]))*arb(lg[i]).exp() for i in range(2)]
x = {-1: Sv[1]*d**(al + 1), 0: Sv[0]*d**al}
g = 1/(d*d*D)
P = {-1: acb(0), 0: acb(1)}; Q = {-1: acb(1), 0: acb(0)}
for j in range(0, JJ):
    beta = al - j
    A, B = BC.coef(beta, d, g, arb(s))
    x[j + 1] = A*x[j] + B*x[j - 1]; P[j + 1] = A*P[j] + B*P[j - 1]; Q[j + 1] = A*Q[j] + B*Q[j - 1]
M0 = float(abs(x[0]).mid()); M1 = float(abs(x[-1]).mid())
df = float(d.mid())
chB = BC.chain_bounds_multi(2.5, nu, nu + 1e-6, M0, M1, JJ, df*(1 - 1e-6), df*(1 + 1e-6), a, a + 1e-6)
for j in (0, 2, 5, 10, 20, 30, 40):
    print(f'j={j}: exact {float(abs(x[j]).mid()):.3g}  sol-frame {float((abs(P[j])*M0 + abs(Q[j])*M1).mid()):.3g}  LemmaB {chB[j+1]:.3g}')
