"""
fullpos/fp_pivot.py (FULLPOS_CERT.md §4.5): size of the §29 signed-remainder terms, with certified moments.
For (J, b, Delta) it evaluates, in Arb,
  S_abs = (w^3/T) sum_j U_j [ |M_j| + (w u_j/(1-u_j/T)) |F_j| ]   (the triangle-inequality majorant),
  M_j = D_j^2 (2Y_j - u_{j-1} L_j),  F_j = sum_{k>j} P_k r^{k-j-1}  (truncated at k <= K: a lower bound on |F_j|),
and the actual H_ang^(2) - H_LR^(2) (k = 2 parts), w = 1 - r.
"""
import sys, os
sys.dont_write_bytecode = True
from flint import arb, fmpq
from math import comb
import fp_common as FC
import fp_poly as FP

Dg = int(sys.argv[1]); K = int(sys.argv[2])
Z, Y = FP.load_moments(Dg, K)
P = [comb(k + 3, 3)**2*(Z[0]/3 if k == 0 else (Z[k] - arb(k)/(k + 3)*Z[k - 1])/(k + 3)) for k in range(K + 1)]
Act = FP.Action(Dg, Z, Y, [fmpq(0)]*(K + 1))
out = []
for (J, b, Delta) in ((30, 0, 1030), (100, 0, 1100), (200, 0, 1200), (100, 0, 3000)):
    T = (Delta - 1)*(Delta - 2); C = J*(J + 1); n = J - b; Ta = arb(T); r = arb(C)/Ta; w = 1 - r
    u = lambda j: (b + j)*(b + j + 1)
    U = arb(1); Sabs = arb(0)
    for j in range(n + 1):
        Lj = Z[0]/3 if j == 0 else (Z[j] - arb(j)/(j + 3)*Z[j - 1])/(j + 3)
        Mj = comb(j + 3, 3)**2*(2*Y[j] - u(j - 1)*Lj)
        Fj = sum((P[k]*r**(k - j - 1) for k in range(j + 1, K + 1)), arb(0))
        Sabs += U*(abs(Mj) + w*u(j)/(1 - arb(u(j))/Ta)*abs(Fj))
        U = U*(C - u(j))/(Ta - u(j))
    Sabs = w**3/Ta*Sabs
    HLR = w**3*sum((P[k]*r**k for k in range(K + 1)), arb(0))
    Hang = Act.H_ang(J, b, T)
    line = (f'J={J} b={b} Delta={Delta} r={float(r.mid()):.4f}: H_ang^(2) = {Hang.str(6)}, H_LR^(2)(trunc K) = {HLR.str(6)}, '
            f'triangle majorant of the remainder >= {Sabs.str(3)}')
    print(line, flush=True); out.append(line)
open(os.path.join(FC.CACHE, f'fp_pivot_{Dg}.txt'), 'w').write('\n'.join(out) + '\n')
