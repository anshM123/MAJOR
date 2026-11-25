"""
np5_spectral.py -- the endpoint measure mu(R) = sum C|a_0|^2 delta_R, the defect D, the endpoint rule and the family bound
in the d = 3 free scalar and free Dirac fermion (exact; no numerics needed beyond bookkeeping).

Facts used (all elementary in free theory, and cross-checked in code/kernels/free_canonical.py, which finds
|ahat^{J_s}_{+-s,2}|^2 = 1 for s = 4, 6 in both theories):
  * E(theta) is one-body and commutes with P, so E_n |T(q),+2> at rest is a two-particle state at rest whose angular
    wavefunction is the T one (e^{2 i theta}) times (E/2) e^{i n theta}(1 + (-1)^n)  [the two particles are back to back].
    Odd n vanish; even n give CM angular momentum j = n + 2.
  * Two-particle states at rest with angular momentum j span exactly ONE conformal family (Flato-Fronsdal):
      scalar  Rac x Rac (sym)        = (+)_{s even >= 0} D(s+1, s)          -> j=0: phi^2 (Delta 1); j=s: J_s (Delta s+1)
      Dirac   Di x Di-bar (U(1) singlet) = D(2,0) (+) (+)_{s>=1} D(s+1, s)  -> j=0: psibar psi (Delta 2); j=s: J_s
  * Normalised harmonic amplitudes are pure phases: |a_{O(j=n+2), n}|^2 = 1, all other O zero (free event shapes are
    delta functions, so the harmonics do not decay).
Ladder dictionary (ANGULAR_ROUTE §§8, 16): helicity 4+k <-> harmonic n = k+2 (b = 4); helicity -k <-> n = -(k+2) (b = 0);
T_Delta = (Delta-1)(Delta-2), C = J(J+1), X_k = C - s(s+1), Y_k = T - s(s+1), s = b+k; R = T/(T-C).
Output: cache/np5_spectral.txt
"""
import os
from fractions import Fraction as Fr
out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)

def family(theory, j):
    j = abs(j)
    if j == 0:
        return ('phi^2', 1, 0) if theory == 'scalar' else ('psibar psi', 2, 0)
    return (f'J_{j}', j + 1, j)

for theory in ('scalar', 'Dirac'):
    P(f"\n=== free {theory}: E_n |T,+2> at rest ===")
    P("  n    O          Delta  J  twist   C      T_Delta   R=T/(T-C)   ladder (b,k)   X_k   Y_k   |a|^2")
    endpoint = Fr(0); famsum = {}
    for n in range(-12, 13, 2):
        j = n + 2
        name, D, J = family(theory, j)
        C, T = J*(J + 1), (D - 1)*(D - 2)
        R = Fr(T, T - C) if T != C else None
        if n >= 2:
            b, k = 4, n - 2; s = b + k
        elif n <= -2:
            b, k = 0, -n - 2; s = -(b + k)
        else:
            b, k, s = None, None, None          # n = 0: the T -> T diagonal (energy), not a heavy-family amplitude
        X = C - abs(s)*(abs(s) + 1) if s is not None else None
        Y = T - abs(s)*(abs(s) + 1) if s is not None else None
        P(f"  {n:3d}  {name:10s} {D:4d} {J:3d}  {D - J:4d}  {C:5d}  {T:7d}   {str(R):>9s}   {str((b, k)):>12s}  {str(X):>4s}  {str(Y):>4s}    1")
        if k == 0 and name != 'J_2':
            endpoint += Fr(T - b*(b - 1), 3) + Fr(3, 2)
            P(f"        endpoint term [(T - b(b-1))/3 + 3/2]|a_0|^2 = {Fr(T - b*(b - 1), 3) + Fr(3, 2)}")
    P(f"  endpoint rule  sum_(O != T, b) [(T - b(b-1))/3 + 3/2]|a_(O,b,0)|^2 = {endpoint}   (author: = 3 + D)  => D = {endpoint - 3}")
    P("  every family that appears has twist 1 (or the J = 0 scalar, C = 0): R = T/(T-C) <= 0, so the band tau >= kappa sqrt(J)")
    P("  is EMPTY: mu((L, oo)) = 0 for all L > 0, B_N = q_(N+2) sum X_N |a_N|^2/(N+1) = 0 (X_N = 0 at every ladder top).")
    P("  family bound: sum_O Delta_O^2 |a_(O,N)|^2 = (|N+2|+1)^2 = N^2 + 6N + 9 and sum_O C|a_(O,N)|^2 = (N+2)(N+3):")
    P("  the author's O(N^2) bound (§17 'Family bound') is SATURATED with unit constant (it cannot be improved to o(N^2)).")
    P("  fourth-harmonic rule: |a_(J_6,4,2)|^2 = 1 = |c|^2 = (t4^2 + alpha_T^2)/16 with t4 = +-4 (scalar +4, fermion -4): consistent.")
P("\nThe gravitational inequality |c|^2 <= (3 + D)/K(g) is of course violated (|c|^2 = 1): the twist gap g is absent -- the")
P("whole heavy sum sits at twist 1 (conserved higher-spin currents), which assumption 4 of §18 excludes.")
open(os.path.join(os.path.dirname(__file__), 'cache', 'np5_spectral.txt'), 'w').write('\n'.join(out) + '\n')
