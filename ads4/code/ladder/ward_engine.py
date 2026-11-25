"""LADDER_ALLJ audit engine: the momentum-space special-conformal action on the bra state <O(p),zeta|, reduced to the rest frame
by Lorentz covariance, applied to a GENERIC zeroth-moment function
      F(p0; zetabar = zeta(A,B), zeta_T = zeta(a,b), z) = sum_m f_m A^{J-m} B^{J+m} a^4 phi_{2-m}(z)     (source helicity +2),
phi_k(z) = homogeneous (degree -2) extension of e^{ik theta}.  Derived here from first principles (LADDER_ALLJ §2):
   K_b = 2(3-Delta) b.d + 2 (p.d)(b.d) - (b.p) d.d  - 2(b.zeta) d_p.d_zeta + 2 (zeta.d_p)(b.d_zeta)        (d = d/dp^mu)
   rest-frame reduction: dF/dp^0 = deg F, dF/dp^i = -A_i F, d_i d_j F = A_i A_j F - delta_ij deg F, d_0 d_i F = -(deg-1) A_i F,
   d_0^2 F = deg(deg-1) F, deg = Delta+1, A_i = total boost generator on (zetabar, zeta_T, z).
Checks:  (1) -K_z F (z = (1,1,0)) equals (W - 18) f  [authors' §7 operator, origin-centred]  -> fixes conventions;
         (2) conservation D_z^mu K_mu F = 0 (Todorov operator, spin -2) -> linear relations on f_m; compared with the §8 recurrence.
Usage: python3 -B ward_engine.py J [J ...] -> cache/ward_engine.txt"""
import sys
import sympy as sp
from common import write

A, B, a, b, z0, z1, z2, D = sp.symbols('A B a b z0 z1 z2 Delta')
ETA = [-1, 1, 1]
Z = [z0, z1, z2]

def spin_ops(x, y):
    Bx = lambda f: sp.Rational(1, 2)*(y*sp.diff(f, x) + x*sp.diff(f, y))
    By = lambda f: sp.I/2*(-y*sp.diff(f, x) + x*sp.diff(f, y))
    Rt = lambda f: -sp.I/2*(x*sp.diff(f, x) - y*sp.diff(f, y))
    return Bx, By, Rt
BxO, ByO, RO = spin_ops(A, B)
BxT, ByT, RT = spin_ops(a, b)
Kxz = lambda f: z1*sp.diff(f, z0) + z0*sp.diff(f, z1)
Kyz = lambda f: z2*sp.diff(f, z0) + z0*sp.diff(f, z2)
Ax = lambda f: sp.expand(BxO(f) + BxT(f) + Kxz(f))
Ay = lambda f: sp.expand(ByO(f) + ByT(f) + Kyz(f))
GEN = {(0, 1): 'x', (1, 0): 'x', (0, 2): 'y', (2, 0): 'y', (1, 2): 'r', (2, 1): 'r'}

def mat_to_spinor(M):
    """decompose an so(2,1) matrix M^rho_sigma into alpha Kx + beta Ky + gamma Rot -> spinor operator on (A,B)"""
    al, be, ga = M[0][1], M[0][2], M[2][1]
    assert M[1][0] == al and M[2][0] == be and M[1][2] == -ga and all(M[i][i] == 0 for i in range(3))
    return lambda f: sp.expand(al*BxO(f) + be*ByO(f) + ga*RO(f))

def Y_op(bvec, mu, sgn=1):
    """vector field  dzeta^rho = sgn*(-2 (b.zeta) eta^{mu rho} + 2 zeta^mu b^rho)  on the O bra polarisation"""
    blow = [ETA[l]*bvec[l] for l in range(3)]
    M = [[sgn*(-2*(ETA[mu] if rho == mu else 0)*blow[sig] + 2*(1 if sig == mu else 0)*bvec[rho]) for sig in range(3)] for rho in range(3)]
    return mat_to_spinor(M)

def Kb(F, bvec, spin_sign=1):
    deg = D + 1
    Ai = {1: Ax, 2: Ay}
    d = {0: sp.expand(deg*F)}
    for i in (1, 2):
        d[i] = sp.expand(-Ai[i](F))
    dd = {(0, 0): sp.expand(deg*(deg - 1)*F)}
    for i in (1, 2):
        dd[(0, i)] = dd[(i, 0)] = sp.expand(-(deg - 1)*Ai[i](F))
        for j in (1, 2):
            dd[(i, j)] = sp.expand(Ai[i](Ai[j](F)) - (deg*F if i == j else 0))
    bp = -bvec[0]
    box = -dd[(0, 0)] + dd[(1, 1)] + dd[(2, 2)]
    out = 2*(3 - D)*sum(bvec[m]*d[m] for m in range(3)) + 2*sum(bvec[m]*dd[(0, m)] for m in range(3)) - bp*box
    out += sum(Y_op(bvec, mu, spin_sign)(d[mu]) for mu in range(3))
    return sp.expand(out)

def phi(k):
    return (z1 + sp.I*z2)**k*z0**(-2 - k) if k >= 0 else (z1 - sp.I*z2)**(-k)*z0**(-2 + k)

def todorov(G, J_spin=-2):
    """D^mu G (upper index), h = 3/2; G homogeneous of degree J_spin in z"""
    h = sp.Rational(3, 2)
    dz = [sp.diff(G, zz) for zz in Z]
    up = [ETA[m]*dz[m] for m in range(3)]
    lap = sum(ETA[m]*sp.diff(G, Z[m], 2) for m in range(3))
    euler = lambda f: sum(Z[m]*sp.diff(f, Z[m]) for m in range(3))
    return [sp.expand((h - 1)*up[m] + euler(up[m]) - sp.Rational(1, 2)*Z[m]*lap) for m in range(3)]

def onshell(expr):
    return sp.expand(expr.subs({a: 1, b: 0}).subs({z0: 1, z1: 1, z2: 0}))

def coeffs(expr, J):
    P = sp.Poly(expr, A, B)
    return {m: sp.factor(P.coeff_monomial(A**(J - m)*B**(J + m))) for m in range(-J, J + 1)}

def W_minus_18(f, J):
    T = (D - 1)*(D - 2); C = J*(J + 1)
    g = lambda m: f.get(m, 0)
    return {m: sp.expand((T + C - 2*m*m - 18)*g(m) - (D - m - 2)*(J + m + 1)*g(m + 1) - (D + m - 2)*(J - m + 1)*g(m - 1))
            for m in range(-J, J + 1)}

def run(J, out):
    pr = lambda x: (print(x, flush=True), out.append(str(x)))
    f = {m: sp.Symbol(f'f{m}'.replace('-', 'm')) for m in range(-J, J + 1)}
    F = sum(f[m]*A**(J - m)*B**(J + m)*a**4*phi(2 - m) for m in range(-J, J + 1))
    target = W_minus_18(f, J)
    for ss in (1, -1):
        Kz = coeffs(onshell(Kb(F, [1, 1, 0], ss)), J)
        ok = all(sp.expand(-Kz[m] - target[m]) == 0 for m in range(-J, J + 1))
        pr(f'J={J}: spin-term sign {ss:+d}:  -K_z F == (W-18) f  (all m, symbolic Delta, generic f): {ok}')
        if ok:
            good = ss
    Kmu = [Kb(F, [1 if n == mu else 0 for n in range(3)], good) for mu in range(3)]
    Dmu = [todorov(Kmu[mu]) for mu in range(3)]
    cons = coeffs(onshell(sum(Dmu[mu][mu] for mu in range(3))), J)     # D^mu K_mu  (D^mu upper, K_mu lower)
    T = (D - 1)*(D - 2); C = J*(J + 1)
    g = lambda m: f.get(m, 0)
    rec = {m: sp.expand(m*(m - 4)*(T + C - 2*m*m)*g(m) - (m - 3)*(m - 4)*(D - m - 2)*(J + m + 1)*g(m + 1)
                        - m*(m - 1)*(D + m - 2)*(J - m + 1)*g(m - 1)) for m in range(-J, J + 1)}
    for m in range(-J, J + 1):
        c = cons[m]
        # which recurrence row is it proportional to?
        match = None
        for mm in range(-J, J + 1):
            if rec[mm] == 0:
                continue
            r = sp.cancel(c/rec[mm])
            if not r.free_symbols - {D}:
                match = (mm, sp.factor(r)); break
        pr(f'   conservation row (coefficient of e_{m}): ' + ('0' if c == 0 else
           (f'= {match[1]} x [authors\' recurrence row m={match[0]}]' if match else f'NOT proportional to a recurrence row: {c}')))
    return out

if __name__ == '__main__':
    out = []
    for J in [int(x) for x in sys.argv[1:]] or [2, 4, 5]:
        run(J, out)
    write('ward_engine.txt', out)
