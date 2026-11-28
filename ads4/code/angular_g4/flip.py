"""Helicity-flip channel of the mixed commutator  <T,-2|[E(th1),E1(th2)]|T,+2>.

Derivation (README §F): same LADDER_ALLJ §4 identity with bra amplitudes u_m = <O,m|E|T,-2> = eta_O v_{-m}
(mirror relation, KERNELS_D3 §5.1, LADDER_ALLJ §1), G' = (i/6) A u, and the same D operator.  Writing the correlator as
e^{2i(th1+th2)} f(phi), phi = th2 - th1, and  i fhat(n) = (1/3) v_{-n}(Av)_n
      + 1/2 [ S_n + (n+3)/2 S_{n-1} - (n+1)/2 S_{n+1} ],   S_m = P_m + P_{-m},  P_m = conj(u_m) v_m.
Contact (KZ: [E,E1] = i delta E):  i fhat(n) = c  for all n.  T family: P_{+-2} = c, A T = 0.
Heavy per-unit-xi weights (xi = eta * alpha * beta, parity-even even-J>=4 states with both ladders), chi_p = sqrt(r0_p r4_{p-4}):
   f_3   = -(1/3) sqrt(r0_3 (C-12)(T-12)) - 2 chi_4                        (rule: sum xi f_3 = -2c)
   f_p   = chi_p (1 + d4_{p-4}/3) + (p+3)/2 chi_{p-1} - (p+1)/2 chi_{p+1}    (p>=4; rule = c)
   f_-p  = chi_p (1 + d0_p/3)     + (3-p)/2 chi_{p+1} + (p-1)/2 chi_{p-1}    (p>=4; rule = c)
b=0 heads (any b=0 state, eta rho0):  F_0: sum eta rho0 (1 + d0_0/3) = c ;  F_1: sum eta rho0 (2 - 4C/3) = 2c.
"""
import numpy as np
from core import ladder, dvals, W
from vec import rmat, dmat, Wv


def chi_mat(J, D, P):
    """chi_p for p = 0..P (p<4 -> 0), vectorised; plus s3 = sqrt(r0_3 (C-12)(T-12))."""
    J = np.asarray(J, float); D = np.asarray(D, float)
    r0 = rmat(np.zeros_like(J), J, D, P + 1)
    r4 = rmat(4 * np.ones_like(J), J, D, P + 1)
    ch = np.zeros((len(J), P + 2))
    for p in range(4, P + 2):
        ch[:, p] = np.sqrt(np.clip(r0[:, p] * r4[:, p - 4], 0, None))
    C = J * (J + 1); T = (D - 1) * (D - 2)
    s3 = np.sqrt(np.clip(r0[:, 3] * (C - 12) * (T - 12), 0, None))
    return ch, s3


def flip_weights(J, D, N):
    """columns: n = 3, 4..N, -4..-N  -> array (S, 1 + 2(N-3)); and the RHS coefficients (c-multiples)."""
    J = np.asarray(J, float); D = np.asarray(D, float)
    ch, s3 = chi_mat(J, D, N + 1)
    d0 = dmat(np.zeros_like(J), J, D, N + 1)
    d4 = dmat(4 * np.ones_like(J), J, D, N + 1)
    cols = [-(1 / 3) * s3 - 2 * ch[:, 4]]; rhs = [-2.0]; labels = [3]
    for p in range(4, N + 1):
        cols.append(ch[:, p] * (1 + d4[:, p - 4] / 3) + (p + 3) / 2 * ch[:, p - 1] - (p + 1) / 2 * ch[:, p + 1]); rhs.append(1.0); labels.append(p)
    for p in range(4, N + 1):
        cols.append(ch[:, p] * (1 + d0[:, p] / 3) + (3 - p) / 2 * ch[:, p + 1] + (p - 1) / 2 * ch[:, p - 1]); rhs.append(1.0); labels.append(-p)
    return np.array(cols).T, np.array(rhs), labels


if __name__ == '__main__':
    # free-field check of the flip rules: currents (chi only at p=J, value eps=1), phi^2 X_0=1, c=1.
    Nn = 30
    X = np.zeros(Nn + 3); X[0] = 1.0
    for Jc in range(4, Nn + 3, 2): X[Jc] = 1.0
    c = 1.0
    S = lambda m: 2 * X[abs(m)] if abs(m) < len(X) else 0.0
    def Tpart(n):
        ST = lambda m: 2 * c if abs(m) == 2 else 0.0
        return 0.5 * (ST(n) + (n + 3) / 2 * ST(n - 1) - (n + 1) / 2 * ST(n + 1))
    out = []
    for n in range(-Nn + 2, Nn - 2):
        heavy = 0.5 * (S(n) + (n + 3) / 2 * S(n - 1) - (n + 1) / 2 * S(n + 1))
        out.append(round(heavy + Tpart(n), 12))
    print('free scalar, i fhat(n) for n=-28..27 (should all be c=1):', set(out))
