"""t5: mode map of the project's C^-_3 kernel into the CLS-D bare-commutator family (LIGHT_INPUTS.md Sec. 3.1, referee M1).

CLS-D light-transform kinematics in the CHMRS frame (SPINNING_E2E.md 3.3, re-derived in LIGHT_INPUTS.md 3.1):
x1 = -y, x3 = y, common point x0 = |y| nhat (nhat _|_ y) => z1 ~ yhat + nhat, z3 ~ yhat - nhat, i.e. the (1,3) polarisations are
ANTIPODAL on the frame circle, theta3 = theta1 + pi.  Hence a bare-commutator functional sees the (1,3) tensor only through
G(phi, phi+pi) = sum_q G_{q1 q3} (-1)^{q3} e^{i(q1+q3)phi}: for total charge M even the weight (-1)^{q3} is symmetric under
q1 <-> q3, for M odd antisymmetric.
The project's C^- uses, at each pair, the parity-odd kernels D5 = eps(z1,z3,yhat) (currents) and A = D5 D1 (stress tensors)
in the preserved-preserved harmonic (M = 0).  We check: D5 and A vanish identically at antipodal polarisations, D1, D1^2 do not;
and the M = 0 modes of D5, A are antisymmetric under q1 <-> q3, hence orthogonal to every antipodal (M even) weight.
"""
import numpy as np
eta = np.diag([-1.0, 1, 1])
yh = np.array([1.0, 0, 0]); ep = np.array([0, 1.0, 0]); eq = np.array([0, 0, 1.0])
z = lambda t: yh + np.cos(t)*ep + np.sin(t)*eq
dot = lambda a, b: a @ eta @ b
eps3 = lambda a, b, c: np.linalg.det(np.array([a, b, c]))
D1 = lambda t1, t3: dot(z(t1), z(t3)) + dot(z(t1), yh)*dot(z(t3), yh)
D5 = lambda t1, t3: eps3(z(t1), z(t3), yh)
kern = {"D1": D1, "D5": D5, "D1^2": lambda a, b: D1(a, b)**2, "A=D5*D1": lambda a, b: D5(a, b)*D1(a, b)}
N = 16
th = 2*np.pi*np.arange(N)/N
print("t5: kernels at antipodal polarisations (CLS-D light-transform kinematics), max over phi:")
for k, f in kern.items():
    anti = max(abs(f(p, p + np.pi)) for p in th)
    G = np.array([[f(a, b) for b in th] for a in th])
    F = np.fft.fft2(G)/N**2                          # F[q1,q3] = coefficient of e^{-i(q1 th1 + q3 th3)}
    modes = {(q1 if q1 <= N//2 else q1-N, q3 if q3 <= N//2 else q3-N): F[q1, q3]
             for q1 in range(N) for q3 in range(N) if abs(F[q1, q3]) > 1e-12}
    M0 = {q: v for q, v in modes.items() if q[0] + q[1] == 0}
    sym = max([abs(v - modes.get((q[1], q[0]), 0)) for q, v in M0.items()] + [0])
    asym = max([abs(v + modes.get((q[1], q[0]), 0)) for q, v in M0.items()] + [0])
    print(f"   {k:8s}: max|K(phi,phi+pi)| = {anti:.2e};  M=0 modes {sorted(M0)}; "
          f"{'symmetric' if sym < 1e-12 else 'ANTISYMMETRIC' if asym < 1e-12 else 'mixed'} under q1<->q3")
print("=> the C^- kernels (D5, A) are invisible to every bare light-ray commutator functional at leading Regge order:")
print("   the project's C^-_3 is NOT the leading moment of a CLS-D bare-commutator (k=3, j'=1) functional.")
