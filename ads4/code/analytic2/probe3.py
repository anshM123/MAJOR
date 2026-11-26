"""probe3: exact S^(alpha - j) via n-recurrence for alpha, alpha+1 (complex128, Re > 0, stable) then the exact beta-recurrence downward
(acb, high precision); true K vs termwise; check one j against the acb n-recurrence."""
import sys, math
import numpy as np
sys.path.insert(0, '.')
from common2 import n_of_L, S_rec
from fastS import S_vec, d_of
from flint import arb, acb, ctx
ctx.prec = 300
L = float(sys.argv[1]); nu = float(sys.argv[2]); a = float(sys.argv[3]); b = int(sys.argv[4]); JJ = int(sys.argv[5])
chk = int(sys.argv[6]) if len(sys.argv) > 6 else -1
s = 2*b + 1
n = n_of_L(L, a, b); 
N = arb(n) + arb(a); D = N*(N + s); d = 1 - arb(n*(n + s))/D
I = acb(0, 1)
al = acb(2.5, nu)
S0, lg = S_vec(n, a, s, [complex(2.5, nu), complex(3.5, nu), complex(4.5, nu)])
Sv = [acb(complex(S0[i]))*arb(lg[i]).exp() for i in range(3)]
# beta recurrence: [D - beta(beta-s)] S(beta-1) = [dD + beta(s-2beta-1)] S(beta) + beta(beta+1) S(beta+1)
S = {0: Sv[0], -1: Sv[1]}     # key -j: S(alpha - j); key -1 means alpha + 1
Sd = {}
Sd[1] = Sv[1]; Sd[0] = Sv[0]
for j in range(0, JJ):
    beta = al - j
    Sd[-(j + 1)] = ((d*D + beta*(s - 2*beta - 1))*Sd[-j] + beta*(beta + 1)*Sd[-j + 1])/(D - beta*(beta - s))
Sh = [Sd[-j]*d**(al - j) for j in range(JJ + 1)]
c = acb(1); a0 = acb(0.5) - I*nu; b0 = acb(1) - 2*I*nu
B = (2*I*nu).gamma()*((acb(0.5) + I*nu).rgamma())**2
tot = acb(0); absum = arb(0)
print(f'n={n} L={float(-d.log()):.4f} nu={nu} a={a} b={b}  kappa={nu/a:.4f}')
for j in range(JJ):
    term = c*Sd[-j]
    tot += term
    absum += abs(c)*abs(Sd[-j])
    if j % 3 == 0 or j < 6:
        print(f' j={j:3d} |S^_j|={float(abs(Sh[j]).mid()):.4e}  ratio={float((abs(Sh[j+1])/abs(Sh[j])).mid()):.4f}  |c_j|d^j={float((abs(c)*d**j).mid()):.3e}  |K partial|={float(abs(2*(B*tot).real).mid()):.4e}')
    c = c*(a0 + j)*(a0 + j)/((b0 + j)*(j + 1))
print('K_true =', float((2*(B*tot).real).mid()), ' |B d^-alpha|=', float(abs(B*d**(-al)).mid()))
print('termwise |B| sum|c_j||S_j| =', float((abs(B)*absum).mid()))
if chk >= 0:
    ctx.prec = 300 + 12*chk
    Sx = S_rec(n, a, s, al - chk, ctx.prec)
    print('check j =', chk, ': n-rec', Sx, ' beta-rec', Sd[-chk])
