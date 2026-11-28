"""Core kinematics for angular_g4.

Conventions follow supplement Sec. sec:angular, ANGULAR_ROUTE §§8-13, ANGULAR_MATCH §§2-3 and code/angmatch/am3_ward.py.
  T = (Delta-1)(Delta-2), C = J(J+1).  Ladders b=0 (helicities m=0,-1,-2,..) and b=4 (m=4,5,...), index j, s=b+j.
  r_j = |a_j|^2/|a_0|^2 = D_j^2 prod_{i<j} (C-s_i(s_i+1))/(T-s_i(s_i+1)),  D_j=(j+1)(j+2)(j+3)/6.
  d_j = (A a)_j / a_j on the ladder = 3[T-(s-1)s]/(j+3) - 3[C-s(s+1)]/(j+1)   (ANGULAR_ROUTE §12; verified in v_dj below).
Per-harmonic mixed-commutator coefficients (exponential Fourier coefficient C^(n), n>=0, of the state-summed
<T,+2|[E(th1),E1(th2)]|T,+2>; am3_ward.py part (a)): a signed harmonic kappa = +-(2+j) with weight w=|a_j|^2, first
moment a_j(Aa)_j = d_j w contributes
     C^(j+2) += (w + d_j w/3)/2,   C^(j+3) += (2 + j/2) w/2,   C^(j+1) += (-j/2) w/2      (same for both ladders).
So per unit |a_0|^2:  h_n = r_{n-2}(1+d_{n-2}/3)/2 + (n+1)/4 r_{n-3} - (n-1)/4 r_{n-1}.
T family: C^(0)=C^(1)=1; the flip entry |c|^2 (harmonic 4, d_T=0) adds gamma_3=-1/2, gamma_4=+1/2, gamma_5=+3/2.
Contact c0 delta only (ANGULAR_MATCH §3: c2=0 by the Ward step)  =>  C^(n) = 1 for every n.
  (M_n)  sum_heavy h_n rho + gamma_n |c|^2 = 1,     n >= 2.
"""
import numpy as np
from fractions import Fraction as Fr

GAMMA = {3: -0.5, 4: 0.5, 5: 1.5}
GAMMA_FR = {3: Fr(-1, 2), 4: Fr(1, 2), 5: Fr(3, 2)}


def allowed(b, J):
    if b == 4:
        return J >= 4
    return J % 2 == 0


def ladder(b, J, D, n, exact=False):
    """r_0..r_n."""
    if exact:
        D = Fr(D); C = Fr(J * (J + 1)); T = (D - 1) * (D - 2); one = Fr(1)
    else:
        C = J * (J + 1.0); T = (D - 1.0) * (D - 2.0); one = 1.0
    r = [one]
    for k in range(n):
        s = b + k
        if k + 1 > J - b:
            r.append(0 * one); continue
        r.append(r[-1] * (one * (k + 4) ** 2 / (k + 1) ** 2) * (C - s * (s + 1)) / (T - s * (s + 1)))
    return r if exact else np.array(r)


def dvals(b, J, D, n, exact=False):
    if exact:
        D = Fr(D); C = Fr(J * (J + 1)); T = (D - 1) * (D - 2)
        return [3 * (T - (b + j - 1) * (b + j)) / (j + 3) - 3 * (C - (b + j) * (b + j + 1)) / Fr(j + 1) for j in range(n + 1)]
    C = J * (J + 1.0); T = (D - 1.0) * (D - 2.0)
    j = np.arange(n + 1); s = b + j
    return 3 * (T - (s - 1) * s) / (j + 3) - 3 * (C - s * (s + 1)) / (j + 1)


def W(b, D, exact=False):
    if exact:
        D = Fr(D); return ((D - 1) * (D - 2) - b * (b - 1)) / 3 + Fr(3, 2)
    return ((D - 1.0) * (D - 2.0) - b * (b - 1)) / 3.0 + 1.5


def hvec(b, J, D, N, exact=False):
    """h_n for n = 0..N (h_0 = h_1 = 0 for heavy states)."""
    r = ladder(b, J, D, N + 1, exact)
    d = dvals(b, J, D, N + 1, exact)
    z = Fr(0) if exact else 0.0
    h = [z] * (N + 1)
    for n in range(2, N + 1):
        v = r[n - 2] * (1 + d[n - 2] / 3) / 2 + (Fr(n + 1, 4) if exact else (n + 1) / 4.0) * (r[n - 3] if n >= 3 else z) \
            - (Fr(n - 1, 4) if exact else (n - 1) / 4.0) * r[n - 1]
        h[n] = v
    return h if exact else np.array(h)


if __name__ == '__main__':
    # v_dj: check d_j against the tridiagonal A on the ladder vector
    rng = np.random.default_rng(0)
    worst = 0
    for t in range(300):
        b = int(rng.choice([0, 4])); J = int(rng.integers(b + 1, 60)); D = J + 2 + 10 ** rng.uniform(-1, 3)
        C = J * (J + 1.0); T = (D - 1) * (D - 2)
        n = J - b + 1
        a = np.sqrt(ladder(b, J, D, n - 1))
        s = np.arange(b, J + 1, dtype=float)
        A = np.diag(T + C - 2 * s ** 2)
        for i in range(n - 1):
            A[i, i + 1] = A[i + 1, i] = -np.sqrt((C - s[i] * (s[i] + 1)) * (T - s[i] * (s[i] + 1)))
        Aa = A @ a
        d = dvals(b, J, D, n - 1)
        m = a > 1e-150
        worst = max(worst, np.max(np.abs(Aa[m] / a[m] - d[m]) / (np.abs(d[m]) + T)))
    print('d_j formula vs tridiagonal A: max rel err', worst)
    # free field: conserved currents J (Delta=J+1, d=0), |a|^2=1 at m=+-J, phi^2 (J=0,Delta=1), |c|^2=1
    Nn = 40
    Ch = np.zeros(Nn + 4)
    def harm(kappa, w, dk=0.0):
        for nn, cf in ((kappa, w + dk / 3), (kappa + 1, (1 + kappa / 2) * w), (kappa - 1, (1 - kappa / 2) * w)):
            if abs(nn) < len(Ch):
                Ch[abs(nn)] += cf if nn == 0 else cf / 2
    harm(0, 1.0); harm(4, 1.0)          # T diag and flip
    harm(2, 1.0)                        # phi^2: m=0 -> kappa=2
    for Jc in range(4, 200, 2):
        harm(2 - Jc, 1.0); harm(2 + Jc, 1.0)
    print('free scalar C^(n), n=0..20:', np.round(Ch[:21], 12))
