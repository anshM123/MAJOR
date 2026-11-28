"""v9: the 3-state family with the J=2 (b=0) state moved BELOW the gap (allowed if it is a light spin-2 single trace,
which the paper's assumptions permit and which A2' does not suppress).  A=(b0,J=2,D_A) free, B=(b4,4,g+4), C=(b4,6,g+6).
Maximise |c|^2 over D_A >= 3 subject to rho_B >= 0.  Compare with D_A = g+2 (no light single traces)."""
from fractions import Fraction as Fr
from math import sqrt
from common import ladder_exact, W

def c2(g, DA):
    A, B, C = (0, 2, DA), (4, 4, Fr(g + 4)), (4, 6, Fr(g + 6))
    rA, rC = ladder_exact(*A, 3), ladder_exact(*C, 3)
    rhoC = Fr(1); rhoA = rhoC * rC[1] / rA[1]; rhoB = rhoA - rhoC
    if rhoB < 0: return None
    s = Fr(3) / (rhoA * W(0, A[2]) + rhoB * W(4, B[2]) + rhoC * W(4, C[2]))
    return s * (rhoC * rC[2] - rhoA * rA[2])

for g in [20, 40, 100, 1000, 10 ** 5]:
    best = (Fr(0), None)
    for i in range(1, 2001):
        DA = Fr(3) + Fr(i, 2000) * (g - 1)
        v = c2(g, DA)
        if v is not None and v > best[0]: best = (v, DA)
    v0 = c2(g, Fr(g + 2))
    print(f"g={g}: gap config tau g^3 = {4*sqrt(v0)*g**3:.2f};  best with light J=2 at D_A={float(best[1]):.2f} "
          f"(D_A/g={float(best[1])/g:.3f}): tau g^3 = {4*sqrt(best[0])*g**3:.2f}  (|c|^2 ratio {float(best[0]/v0):.3f})")
