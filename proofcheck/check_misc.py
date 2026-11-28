"""Misc. independent checks:
1. Renormalised Crofton endpoint integral (SM Prop. S4.4): m{u>u0, v>v0, uv<Y^2} with dm = du dv/(u+v)^2
   equals the closed form and -> log(Y/d_p); Remark S4.5(d) end constant.
2. Explicit invariant strip at a rotating horizon (SM Sec. S10, Eqs. strip, stripmetric, separation identity).
3. Lemma N3 (d_n increase to the closed length) for vacuum data on a hyperbolic cylinder (sanity check of the
   cyclic-order bookkeeping), D(x,y)=2log|x-y|-log|x|-log|y| for h(z)=lam z.
4. 3d Jacobi coefficients (SM S17): rho = 1 - T(k,l), <R(k,t)t,k> = T(k,k) from the 3d Riemann/Schouten form.
"""
import mpmath as mp
import sympy as sp

# 1 ---------------------------------------------------------------------------------------------
mp.mp.dps = 30
for (u0, v0) in [(0.3, 0.7), (1.0, 0.0), (0.0, 2.0), (0.05, 0.05)]:
    d = u0 + v0
    for Y in [10.0, 1e3, 1e5]:
        if v0 > 0:
            num = mp.quad(lambda u: 1 / (u + v0) - u / (u * u + Y * Y), [u0, Y * Y / v0])
            closed = mp.log((Y * Y / v0 + v0) / d) - mp.log((Y ** 4 / v0 ** 2 + Y * Y) / (u0 ** 2 + Y * Y)) / 2
        else:
            num = mp.quad(lambda u: 1 / u - u / (u * u + Y * Y), [u0, mp.inf])
            closed = mp.log(1 + Y * Y / u0 ** 2) / 2
        print(f'u0={u0} v0={v0} Y={Y:g}: numeric={float(num):.12f} closed={float(closed):.12f} '
              f'log(Y/d)={float(mp.log(Y/d)):.12f}')
# end constant: m{a<q-, b>q+, b-a<2Y}
d = 1.3
for Y in [1e2, 1e4]:
    num = mp.quad(lambda s: (s - d) / s ** 2, [d, 2 * Y])
    print('end-constant check:', float(num - (mp.log(2 * Y / d) - 1)))

# 2 ---------------------------------------------------------------------------------------------
s, r, s2, r2, a, b = sp.symbols('s r s2 r2 a b', positive=True)
L = (a + b) / 2
q = (a - b) / (a + b)
c = lambda S: sp.Matrix([[0, sp.exp(S)], [-sp.exp(-S), 0]])
n = lambda S: sp.Matrix([[sp.exp(q * S), 0], [0, -sp.exp(-q * S)]])
X = lambda S, R: sp.cosh(R) * c(S) + sp.sinh(R) * n(S)
print('det X - 1 =', sp.simplify(X(s, r).det() - 1))
hL = sp.diag(sp.exp(a / 2), sp.exp(-a / 2))
hR = sp.diag(sp.exp(b / 2), sp.exp(-b / 2))
print('equivariance residual =', sp.simplify(hL * X(s, r) * hR.inv() - X(s + L, r)))
Xs, Xr = X(s, r).diff(s), X(s, r).diff(r)
Nf = lambda P, Q: (P.adjugate() * Q).trace()
# metric -det(dX): g(V,W) = -(1/2) N(V,W)
gss = sp.simplify(-Nf(Xs, Xs) / 2)
grr = sp.simplify(-Nf(Xr, Xr) / 2)
grs = sp.simplify(-Nf(Xr, Xs) / 2)
print('g_ss - [1+(1-q^2)sinh^2 r] =', sp.simplify(gss - (1 + (1 - q ** 2) * sp.sinh(r) ** 2)),
      '; g_rr =', grr, '; g_rs =', grs)
Dl = sp.symbols('Delta')
sep = Nf(X(s, r), X(s - Dl, r2)) / 2
claim = sp.cosh(r - r2) * sp.cosh(q * Dl) + sp.cosh(r) * sp.cosh(r2) * (sp.cosh(Dl) - sp.cosh(q * Dl))
print('separation identity residual =', sp.simplify((sep - claim).rewrite(sp.exp)))

# 3 ---------------------------------------------------------------------------------------------
lam = 3.0
Dv = lambda x, y: 2 * mp.log(abs(x - y)) - mp.log(abs(x)) - mp.log(abs(y))
x0, y0 = 1.0, -2.0
dn = [Dv(x0, lam ** n * y0) - Dv(x0, lam ** (n - 1) * y0) for n in range(1, 12)]
print('d_n:', [round(float(v), 6) for v in dn], ' log(lam)=', float(mp.log(lam)))
print('monotone:', all(dn[i] <= dn[i + 1] + 1e-15 for i in range(len(dn) - 1)))

# 4 ---------------------------------------------------------------------------------------------
# frame tau (spacelike unit), k, l null, <k,l>=-1. Use components in basis (tau,k,l).
gm = sp.Matrix([[1, 0, 0], [0, 0, -1], [0, -1, 0]])
Tsym = sp.Matrix(3, 3, lambda i, j: sp.Symbol('T%d%d' % (min(i, j), max(i, j))))
gi = gm.inv()
Ttr = sum(gi[i, j] * Tsym[i, j] for i in range(3) for j in range(3))
# Einstein with Lambda=-1, 8piG=1: Ric = T - (tr T + 2) g ; Schouten P = Ric - R/4 g
Ric = Tsym - (Ttr + 2) * gm
Rs = sum(gi[i, j] * Ric[i, j] for i in range(3) for j in range(3))
P = Ric - Rs / 4 * gm
Rm = lambda A, B, C, D: gm[A, C] * P[B, D] + gm[B, D] * P[A, C] - gm[A, D] * P[B, C] - gm[B, C] * P[A, D]
# check against constant curvature -1 when T=0 and against Ricci contraction
zeroT = {sym: 0 for sym in Tsym.free_symbols}
ricc = sp.Matrix(3, 3, lambda B, D: sum(gi[A, C] * Rm(A, B, C, D) for A in range(3) for C in range(3)))
print('Riemann->Ricci consistency:', sp.simplify(ricc - Ric) == sp.zeros(3, 3))
tau, k, l = 0, 1, 2
# <R(X,Y)Y,X> = Rm(X,Y,X,Y)  (sectional-curvature convention);  <R(k,tau)tau,l> = Rm(k,tau,l,tau)
print('rho = <R(k,t)t,l> =', sp.simplify(Rm(k, tau, l, tau)), ' (claim 1 - T_kl =', 1 - Tsym[1, 2], ')')
print('<R(k,t)t,k> =', sp.simplify(Rm(k, tau, k, tau)), ' Ric_kk =', sp.simplify(Ric[k, k]))
print('vacuum sectional K(tau,k-l...) check: Rm(t,k,t,l) at T=0 =', Rm(tau, k, tau, l).subs(zeroT))
