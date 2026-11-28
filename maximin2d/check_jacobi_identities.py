"""Numerical check of the 3D curvature identities used in Prop. J, in Vaidya-AdS3 at a random point
with a random spacelike unit tau and its normal null frame k,l (<k,l>=-1):
  (i)  <R(l,tau)tau,l> = Ric(l,l)          (coupling coefficient = NEC quantity)
  (ii) rho := <R(k,tau)tau,l> = 1 - T(k,l), T := Ric - R/2 g - g  (Einstein, Lambda=-1, 8piG=1)
  (iii) Jacobi system alpha'' = rho alpha + Ric(l,l) beta for J = alpha k + beta l
Convention: <R(X,Y)Z,W> with R(X,Y)Z = nabla_X nabla_Y Z - ...; Jacobi eq J'' = -R(J,tau)tau.
usage: python3 -I check_jacobi_identities.py"""
import sympy as sp, numpy as np
v, z, x = sp.symbols('v z x', real=True)
X = [v, z, x]
m = 1 + sp.Rational(3, 2) * (1 + sp.tanh(v / sp.Rational(3, 10))) / 2
f = 1 - m * z**2
g = sp.Matrix([[-f, -1, 0], [-1, 0, 0], [0, 0, 1]]) / z**2
gi = g.inv()
n = 3
Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(n)) / 2
         for c in range(n)] for b in range(n)] for a in range(n)]
# R^a_{bcd} = d_c Gam^a_{db} - d_d Gam^a_{cb} + Gam^a_{ce}Gam^e_{db} - Gam^a_{de}Gam^e_{cb};  R(X,Y)Z^a = R^a_{bcd} Z^b X^c Y^d
Rud = [[[[sp.diff(Gam[a][d][b], X[c]) - sp.diff(Gam[a][c][b], X[d]) + sum(Gam[a][c][e] * Gam[e][d][b] - Gam[a][d][e] * Gam[e][c][b] for e in range(n))
          for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
pt = {v: 0.13, z: 0.7, x: 0.0}
Rn = np.array([[[[float(Rud[a][b][c][d].subs(pt)) for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)])
gn = np.array(g.subs(pt), dtype=float)
Ric = np.einsum('abad->bd', Rn)                     # Ric_{bd} = R^a_{bad}
Rs = np.einsum('ab,ab', np.linalg.inv(gn), Ric)
T = Ric - 0.5 * Rs * gn - gn
ip = lambda A, B: A @ gn @ B
def Rvec(Xv, Yv, Zv): return np.einsum('abcd,b,c,d->a', Rn, Zv, Xv, Yv)
rng = np.random.default_rng(1)
# random spacelike unit tau, then orthonormal normal frame (n timelike future, N spacelike)
for trial in range(3):
    while True:
        t = rng.normal(size=3)
        if ip(t, t) > 0.1: break
    t /= np.sqrt(ip(t, t))
    # timelike vector: -d_z direction? use Gram-Schmidt on candidates
    cand = [np.array([1.0, -0.5, 0.0]), np.array([0.0, 0.0, 1.0]), np.array([0.3, 1.0, 0.2])]
    basis = [t]
    for c in cand:
        w = c.copy()
        for b in basis: w = w - ip(w, b) / ip(b, b) * b
        if abs(ip(w, w)) > 1e-6: basis.append(w / np.sqrt(abs(ip(w, w))))
    nn = [b for b in basis[1:] if ip(b, b) < 0][0]; NN = [b for b in basis[1:] if ip(b, b) > 0][0]
    k = (nn + NN) / np.sqrt(2); l = (nn - NN) / np.sqrt(2)
    print(f"trial {trial}: <k,l>={ip(k,l):+.3f} <k,k>={ip(k,k):+.1e}")
    print("  (i)  <R(l,t)t,l> =", round(ip(Rvec(l, t, t), l), 10), "  Ric(l,l) =", round(l @ Ric @ l, 10))
    print("       <R(k,t)t,k> =", round(ip(Rvec(k, t, t), k), 10), "  Ric(k,k) =", round(k @ Ric @ k, 10))
    print("  (ii) rho=<R(k,t)t,l> =", round(ip(Rvec(k, t, t), l), 10), "  1-T(k,l) =", round(1 - k @ T @ l, 10))
