"""v8: (b) compare the two ladders on a heavy state.  r_k^(4)/r_k^(0) at fixed J as T_D -> infinity tends to
prod_{j<k} (J-4-j)(J+5+j)/((J-j)(J+1+j)), not 1: the ladders share the power T^-k but not the coefficient unless J >> k.
Also checks the closed Pochhammer form r_k = binom(k+3,3)^2 (b-J)_k (J+b+1)_k / ((b+2-D)_k (D-1+b)_k)
against the paper's step ratio, and the paper's |a_2|^2/|a_0|^2 = 100 (C-20)(C-30)/((T-20)(T-30)) on b=4."""
from fractions import Fraction as Fr
from math import comb
from common import ladder_exact

def poch(x, k):
    p = Fr(1)
    for i in range(k): p *= (x + i)
    return p

for (b, J, D) in [(0, 10, Fr(57, 2)), (4, 12, Fr(41)), (4, 9, Fr(1001, 7))]:
    r = ladder_exact(b, J, D, 8)
    for k in range(9):
        if k > J - b: break
        alt = comb(k + 3, 3) ** 2 * poch(b - J, k) * poch(J + b + 1, k) / (poch(b + 2 - D, k) * poch(D - 1 + b, k))
        assert alt == r[k], (b, J, D, k)
C = lambda J: J * (J + 1); T = lambda D: (D - 1) * (D - 2)
for J, D in [(6, Fr(37)), (10, Fr(53, 3))]:
    assert ladder_exact(4, J, D, 2)[2] == 100 * Fr((C(J) - 20) * (C(J) - 30)) / ((T(D) - 20) * (T(D) - 30))
print("Pochhammer form and paper's k=2 ratio: exact agreement")
for J in [6, 8, 12, 20, 50, 200]:
    D = Fr(10 ** 6)
    r0, r4 = ladder_exact(0, J, D, 3), ladder_exact(4, J, D, 3)
    print(f"  J={J:4d}, Delta=1e6:  r1^(4)/r1^(0) = {float(r4[1]/r0[1]):.4f}   r2^(4)/r2^(0) = {float(r4[2]/r0[2]):.4f}")
