"""v5: print the primal (amplitude) configuration behind a v3 dual optimum (HiGHS marginals), to see which
states saturate the LP; used for the light-column variant.  usage: python3 v5_primal_diag.py g N [light]"""
import sys, numpy as np
import v3_lp_hp as V
from scipy.optimize import linprog
g = int(sys.argv[1]); N = int(sys.argv[2]); light = 'light' in sys.argv
M = V.Mmat(N)
MT = V.arb_mat([[V.arb(M[n][k].numerator) / M[n][k].denominator for n in range(N + 1)] for k in range(N + 1)])
S = V.candidates(g, light)
idx = [i for i, s in enumerate(S) if s[1] <= 3 * g + 10 and float(s[2]) - s[1] <= 30 * g]
S = [S[i] for i in idx]
B = V.basis_rows(S, N, MT)
a_eq = np.array([float(M[n][2]) for n in range(N + 1)])
A = np.hstack([B, -np.ones((len(S), 1))]); c = np.zeros(N + 2); c[-1] = 1
res = linprog(c, A_ub=A, b_ub=np.zeros(len(S)), A_eq=np.append(a_eq, 0)[None, :], b_eq=[100.0],
              bounds=[(None, None)] * (N + 1) + [(0, None)], method='highs')
lam = res.x[-1]; y = -res.ineqlin.marginals          # y_i = rho_i W_i (up to overall normalisation)
y = 3 * y / y.sum()
print(f"g={g} N={N} light={light}: |c|^2 = {3*lam:.6e}")
for i in np.argsort(-y)[:10]:
    if y[i] > 1e-6:
        b, J, D = S[i]
        print(f"   b={b} J={J} Delta={float(D):.4f} twist={float(D)-J:.4f}  share of endpoint budget y={y[i]:.4f}")
