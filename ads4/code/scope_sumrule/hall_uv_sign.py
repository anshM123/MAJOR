"""
Sign of the UV (OPE) constant of the Hall correlator relative to the Hall viscosity, for the parity-odd cubic
coupling Q = C C C~ on the AdS4 black brane, at first order in its coupling lam.  [computed]

Conventions are exactly those of ../scope_thermal/hall/hall_ef.py (ORIENT = +1, eps_{vrxy} = +r^2):
    ds^2 = -f dv^2 + 2 dv dr + r^2 (e^{eps H})_ij dx^i dx^j,  H = [[phi, chi], [chi, -phi]],  f = r^2 - 1/r,
    S = (1/16piG) Int sqrt(-g) (R + 6 + lam Q).
phi sources T1 = (Txx - Tyy)/2, chi sources T2 = Txy (delta S = Int (T1 phi + T2 chi)).
h(w) := G^R_{T2 T1}(w) = <T^xy> response to phi_0 e^{-i w v}.  In the same conventions hall_ef.py finds
    <T^xy> = -p chi - eta chi_v + eta_H phi_v   with 4 pi eta_H/s = +108 lam,
i.e. h(w) = -i w eta_H + O(w^2); the IR identity (1/pi) Int_0^inf (d rho_+ - d rho_-)/w^2 dw = eta_H holds with the
same eta_H (contour argument in README).  The even channel g(w) = P + i eta w + ... has rho/w -> +eta (CDP sign).

We compute, at Euclidean frequency w = i Omega (e^{-iwv} = e^{Omega v}),
    h(i Omega) = Pi_min(oo)/(16 pi G),  Pi_min = -r^2 (f X_r + Omega X),
where chi = lam X(r) e^{Omega v} solves the Einstein shear operator with the source from the Q cross terms acting on
the zeroth-order regular solution phi = Phi(r) e^{Omega v}, Phi(oo) = 1, and X has no source at the boundary.
  small Omega:   h(i Omega) -> eta_H Omega            (check: Pi(oo)/Omega -> 108)
  large Omega:   h(i Omega) -> J_H  (OPE constant),   S_H = J_H / eps = Pi(oo)/2   (eps = 2P = 2/(16 pi G)).
The collaborator's combination Z = V + (12/5)(S_E - 3/8 + i S_H), V = (4 pi eta/s - 1) - i 4 pi eta_H/s, cancels the
cubic imaginary part iff (12/5) S_H = 4 pi eta_H/s, i.e. S_H = +45 lam here; it doubles it if S_H = -45 lam.

Method: Chebyshev collocation in z = 1/r on [0, 1] (regularity at the horizon is automatic; Dirichlet at z = 0).
The Q cross Lagrangian is the one printed by hall_ef.py (out_hall_ef_orient+1.txt), times sqrt(-g) = r^2.
"""
import numpy as np
import sympy as sp

v, r, z = sp.symbols('v r z', positive=True)
Om = sp.Symbol('Omega', positive=True)
names = ['p', 'pv', 'pr', 'pvv', 'pvr', 'prr', 'q', 'qv', 'qr', 'qvv', 'qvr', 'qrr']
Js = {k: sp.Symbol(k) for k in names}
Qcross = ("-54*pr*qv/r**6 + 24*pr*qvr*(r - 1)*(r**2 + r + 1)/r**5 + 24*pr*qvv/r**4 - 18*prr*qv/r**5 "
          "+ 12*prr*qvr*(r - 1)*(r**2 + r + 1)/r**4 + 12*prr*qvv/r**3 + 54*pv*qr/r**6 + 18*pv*qrr/r**5 "
          "- 24*pvr*qr*(r - 1)*(r**2 + r + 1)/r**5 - 12*pvr*qrr*(r - 1)*(r**2 + r + 1)/r**4 - 24*pvv*qr/r**4 "
          "- 12*pvv*qrr/r**3")
L = r**2*sp.sympify(Qcross, locals={**Js, 'r': r})

phi = sp.Function('phi')(v, r)
rev = {Js['p']: phi, Js['pv']: sp.diff(phi, v), Js['pr']: sp.diff(phi, r), Js['pvv']: sp.diff(phi, v, 2),
       Js['pvr']: sp.diff(phi, v, r), Js['prr']: sp.diff(phi, r, 2)}
E = 0
for k, ops in (('q', ()), ('qv', (v,)), ('qr', (r,)), ('qvv', (v, v)), ('qvr', (v, r)), ('qrr', (r, r))):
    dl = sp.diff(L, Js[k]).subs(rev)
    for var in ops:
        dl = -sp.diff(dl, var)
    E += dl
E = sp.expand(E)
# O(w) check against hall_ef: terms phi_v * sigma(r)
PV = sp.Derivative(phi, v)
sig = sum(t/PV for t in sp.Add.make_args(E)
          if [a for a in t.atoms(sp.Derivative) if a.expr == phi] == [PV] and not (t/PV).has(phi))
print('sigma(r) from this EL =', sp.simplify(sig), '  (hall_ef: 432/r**5);  int_1^oo =', sp.integrate(sig, (r, 1, sp.oo)))

# substitute phi = Phi(z) e^{Omega v}, z = 1/r
Phi = sp.Function('Phi')
Phir = Phi(1/r)
Es = E.subs(phi, Phir*sp.exp(Om*v)).doit()
Es = sp.simplify(Es*sp.exp(-Om*v))
Es = Es.subs(r, 1/z)
Es = sp.expand(sp.simplify(Es))
dPhi = [Phi(z)] + [sp.Derivative(Phi(z), (z, k)) for k in range(1, 5)]
P0, P1, P2, P3, P4 = sp.symbols('P0:5')
Es = Es.subs({dPhi[4]: P4, dPhi[3]: P3, dPhi[2]: P2, dPhi[1]: P1}).subs(Phi(z), P0)
Es = sp.expand(Es)
coef = [sp.simplify(Es.coeff(Pk)) for Pk in (P0, P1, P2, P3, P4)]
assert sp.simplify(Es - sum(c*Pk for c, Pk in zip(coef, (P0, P1, P2, P3, P4)))) == 0
print('source S(z) = sum_k c_k(z, Omega) Phi^(k)(z):')
for k, c in enumerate(coef):
    print('  c_%d =' % k, sp.factor(c))
cfun = [sp.lambdify((z, Om), c, 'numpy') for c in coef]


def cheb(N):
    """Chebyshev points on [0,1] (x_0 = 1 ... x_N = 0) and differentiation matrix (Trefethen)."""
    k = np.arange(N + 1)
    x = np.cos(np.pi*k/N)
    c = np.hstack([2, np.ones(N - 1), 2])*(-1)**k
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1/c)/(dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    zz = (x + 1)/2
    return zz, 2*D


def solve(Omv, N=48):
    zz, D = cheb(N)
    D2 = D @ D
    I = np.eye(N + 1)
    # operator z(1-z^3) y'' - (2 + z^3 + 2 Om z) y' + 2 Om y   (= z * [Einstein shear operator in r])
    A = np.diag(zz*(1 - zz**3)) @ D2 - np.diag(2 + zz**3 + 2*Omv*zz) @ D + 2*Omv*I
    ib = N  # index of z = 0
    # Phi: Dirichlet Phi(0) = 1
    A1 = A.copy(); b1 = np.zeros(N + 1)
    A1[ib, :] = 0; A1[ib, ib] = 1; b1[ib] = 1
    Ph = np.linalg.solve(A1, b1)
    ders = [Ph]
    for _ in range(4):
        ders.append(D @ ders[-1])
    S = sum(cfun[k](zz, Omv)*ders[k] for k in range(5))
    # X: z*L[X] = -z*S, X(0) = 0
    b2 = -zz*S
    A2 = A.copy()
    A2[ib, :] = 0; A2[ib, ib] = 1; b2[ib] = 0
    X = np.linalg.solve(A2, b2)
    g = (1 - zz**3)*(D @ X) - Omv*X          # Pi = g / z^2
    g1 = D @ g; g2 = D @ g1
    return g[ib], g1[ib], g2[ib]/2, Ph


if __name__ == '__main__':
    print('\nPi(oo) = g''(0)/2 with g = (1-z^3) X_z - Omega X (requires g(0) = g\'(0) = 0)')
    for Omv in (1e-3, 1e-2):
        for N in (40, 60):
            g0, g1, Pi, _ = solve(Omv, N)
            print('Omega=%-7g N=%d  g(0)=%.1e g\'(0)=%.1e  Pi/Omega = %.8f   (expect 4 pi eta_H/s = 108)' % (Omv, N, g0, g1, Pi/Omv))
    print()
    rows = []
    for Omv in (0.5, 1, 2, 5, 10, 20, 40, 80, 160):
        vals = []
        for N in (60, 90):
            g0, g1, Pi, _ = solve(Omv, N)
            vals.append(Pi)
        print('Omega=%-6g Pi(oo) = %.10f (N=60)  %.10f (N=90)   S_H = Pi/2 = %.6f' % (Omv, vals[0], vals[1], vals[1]/2))
        rows.append((Omv, vals[1]))
    # Richardson in 1/Omega^3 (next OPE term ~ Omega^-3 from Delta = 6)
    (o1, p1), (o2, p2) = rows[-2], rows[-1]
    ext = (p2*o2**3 - p1*o1**3)/(o2**3 - o1**3)
    print('\nlarge-Omega extrapolation: S_H = J_H/eps = %.6f per unit lam   (|alpha_T|/32 predicts |S_H| = 45)' % (ext/2))
    print('4 pi eta_H/s = +108 lam (hall_ef, same conventions).  Ratio S_H / (4 pi eta_H/s) = %.6f' % (ext/2/108))
    print('=> Z = V + (12/5)(S_E - 3/8 + i S_H): imaginary cubic part', 'CANCELS' if ext > 0 else 'ADDS (doubles)')
