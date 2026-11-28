"""Check the flip-channel Fourier formula against the defining LADDER_ALLJ §4 expression with bra amplitudes u
(numerically, random complex u, v and random real symmetric A on helicities -J..J):
  LHS(th1,th2) = sum_m [ conj(G'_m(th1)) F_m(th2) - conj(F'_m(th2)) G_m(th1) - (i/2)( conj(F'_m(th1)) (D F_m)(th2) + (D conj(F'_m))(th2) F_m(th1) ) ]
  F_m = v_m e^{i(2-m)th}, F'_m = u_m e^{i(-2-m)th}, G = (i/6) A F, G' = (i/6) A F',  D = 1 + 2cos(phi) + sin(phi) d/dth2, phi = th2-th1.
Claim: LHS = -i e^{2i(th1+th2)} sum_n X_n e^{i n phi} with (X_n =: i fhat(n), contact => X_n = c for all n)
  i fhat(n) = (1/6)(Q_{-n} + Q'_n) + 1/2[S_n + (n+3)/2 S_{n-1} - (n+1)/2 S_{n+1}],
  P_m = conj(u_m) v_m, S_m = P_m + P_{-m}, Q_m = conj((Au)_m) v_m, Q'_m = conj(u_m)(Av)_m.
(The general-u version keeps P_m and P_{-m} separately; check that form.)"""
import numpy as np
rng = np.random.default_rng(3)
J = 6; ms = np.arange(-J, J + 1); n_ = len(ms)
u = rng.normal(size=n_) + 1j * rng.normal(size=n_); v = rng.normal(size=n_) + 1j * rng.normal(size=n_)
A = rng.normal(size=(n_, n_)); A = A + A.T
Au, Av = A @ u, A @ v
def lhs(t1, t2):
    phi = t2 - t1
    F = lambda t: v * np.exp(1j * (2 - ms) * t)
    Fp = lambda t: u * np.exp(1j * (-2 - ms) * t)
    G = lambda t: 1j / 6 * Av * np.exp(1j * (2 - ms) * t)
    Gp = lambda t: 1j / 6 * Au * np.exp(1j * (-2 - ms) * t)
    k = 2 - ms; kp = 2 + ms     # harmonics of F(th2) and conj(F')(th2)
    DF = F(t2) * (1 + (1 + k / 2) * np.exp(1j * phi) + (1 - k / 2) * np.exp(-1j * phi))
    DFp = np.conj(Fp(t2)) * (1 + (1 + kp / 2) * np.exp(1j * phi) + (1 - kp / 2) * np.exp(-1j * phi))
    return np.sum(np.conj(Gp(t1)) * F(t2) - np.conj(Fp(t2)) * G(t1) - 0.5j * (np.conj(Fp(t1)) * DF + DFp * F(t1)))
P = {m: np.conj(u[i]) * v[i] for i, m in enumerate(ms)}
Q = {m: np.conj(Au[i]) * v[i] for i, m in enumerate(ms)}
Qp = {m: np.conj(u[i]) * Av[i] for i, m in enumerate(ms)}
g = lambda d, m: d.get(m, 0)
def ifhat(n):
    return (1 / 6) * (g(Q, -n) + g(Qp, n)) + 0.5 * (g(P, n) + g(P, -n) + (n + 3) / 2 * (g(P, n - 1) + g(P, 1 - n)) - (n + 1) / 2 * (g(P, n + 1) + g(P, -n - 1)))
err = 0
for t1, t2 in rng.uniform(0, 2 * np.pi, size=(20, 2)):
    phi = t2 - t1
    rhs = -1j * np.exp(2j * (t1 + t2)) * sum(ifhat(n) * np.exp(1j * n * phi) for n in range(-2 * J - 3, 2 * J + 4))
    err = max(err, abs(lhs(t1, t2) - rhs))
print('max |LHS - claimed| over random points:', err)
