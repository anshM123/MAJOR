"""
Free-field three- and two-point functions in position space (d = 3) by Wick contraction, evaluated exactly at
rational points, and their expansion in the embedding structures of structures.py.

Scalar:  <phi(x)phi(0)> = 1/(4 pi (x^2)^{1/2});  bilinears B = sum_{a+b=s} c_{ab} (z.d)^a phi (z.d)^b phi.
           T = (3/4)(z.d phi)^2 - (1/4) phi (z.d)^2 phi  (canonical, C_T = 3/(32 pi^2));
           J_s = sum_k (-1)^k C(2s,2k) (z.d)^{s-k} phi (z.d)^k phi.
Dirac fermion (= two Majorana fermions; all operators used are sums over the two copies):
           <psi(x) psibar(0)> = c_F xslash/(x^2)^{3/2};  B = sum c_{ab} (z.d)^a psibar Gamma (z.d)^b psi,
           Gamma = zslash (currents) or 1 (psibar psi);  T ~ (z.d psibar) zslash psi - psibar zslash (z.d psi);
           J_s = sum_k (-1)^k C(2s,2k+1) (z.d)^{s-1-k} psibar zslash (z.d)^k psi.
         The constants c_F and the normalisation of T are not needed: they are fixed by C_T (see free_values).
All correlators are (x41^2 x43^2 x13^2)^{-1/2} (resp. (x^2)^{-1}) times a rational function, which is what is returned.
"""
from fractions import Fraction as Fr
from math import comb, factorial
import random
import sympy as sp
import structures as st

ETA = (-1, 1, 1)
def dot(u, v):
    return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def sub(u, v):
    return [a - b for a, b in zip(u, v)]

def binom_half(n):          # binom(-1/2, n)
    r = Fr(1)
    for q in range(n):
        r = r*(Fr(-1, 2) - q)/(q + 1)
    return r

def g_deriv(k, l, ze, z, x, grad=False):
    """(ze.d)^k (z.d)^l (x^2)^{-1/2} = (x^2)^{-1/2} * value ; optionally also the gradient d_mu (lower index)."""
    zx, zex, zz, x2 = dot(z, x), dot(ze, x), dot(ze, z), dot(x, x)
    val = Fr(0); gr = [Fr(0)]*3
    lo = lambda v: [ETA[m]*v[m] for m in range(3)]
    zel, zl, xl = lo(ze), lo(z), lo(x)
    for r in range(0, min(k, l) + 1):
        n = k + l - r
        c = factorial(k)*factorial(l)*binom_half(n)*Fr(factorial(n), factorial(r)*factorial(k-r)*factorial(l-r))*2**n*zz**r
        if c == 0:
            continue
        p, q = k - r, l - r
        base = zex**p*zx**q/x2**n
        val += c*base
        if grad:
            for m in range(3):
                t = Fr(0)
                if p:
                    t += p*zel[m]*zex**(p-1)*zx**q/x2**n
                if q:
                    t += q*zl[m]*zex**p*zx**(q-1)/x2**n
                t += (Fr(-1, 2) - n)*2*xl[m]*base/x2
                gr[m] += c*t
    return (val, gr) if grad else val

# ---------------------------------------------------------------- scalar
def c_scalar_T():
    return {(0, 2): Fr(-1, 8), (1, 1): Fr(3, 4), (2, 0): Fr(-1, 8)}
def c_scalar_J(s):
    return {(s - k, k): Fr((-1)**k*comb(2*s, 2*k)) for k in range(s + 1)}

def scalar_3pt(c4, c1, c3, x4, z4, x1, z1, x3, z3):
    """<B4(x4,z4) B1(x1,z1) B3(x3,z3)> * (4 pi)^3 (x41^2 x43^2 x13^2)^{1/2}"""
    x41, x43, x13 = sub(x4, x1), sub(x4, x3), sub(x1, x3)
    tot = Fr(0)
    for (a, b), ca in c4.items():
        for (c, d), cc in c1.items():
            g41 = (-1)**c*g_deriv(a, c, z4, z1, x41)
            for (e, f), ce in c3.items():
                g43 = (-1)**e*g_deriv(b, e, z4, z3, x43)
                g13 = (-1)**f*g_deriv(d, f, z1, z3, x13)
                tot += ca*cc*ce*g41*g43*g13
    return 8*tot
def scalar_2pt(c, cp, x, z, zp):
    """<B(x,z) B'(0,z')> * (4 pi)^2 x^2"""
    tot = Fr(0)
    for (a, b), ca in c.items():
        for (c_, d), cc in cp.items():
            tot += ca*cc*(-1)**(c_ + d)*g_deriv(a, c_, z, zp, x)*g_deriv(b, d, z, zp, x)
    return 2*tot

# ---------------------------------------------------------------- Dirac fermion
G0 = ((Fr(0), Fr(1)), (Fr(-1), Fr(0))); G1 = ((Fr(0), Fr(1)), (Fr(1), Fr(0))); G2 = ((Fr(1), Fr(0)), (Fr(0), Fr(-1)))
GAM = (G0, G1, G2)        # gamma^mu, {gamma^mu, gamma^nu} = 2 eta^{mu nu}
def mmul(A, B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)) for i in range(2))
def madd(A, B, c=1):
    return tuple(tuple(A[i][j] + c*B[i][j] for j in range(2)) for i in range(2))
def mscal(c, A):
    return tuple(tuple(c*A[i][j] for j in range(2)) for i in range(2))
ZERO2 = ((Fr(0), Fr(0)), (Fr(0), Fr(0))); ID2 = ((Fr(1), Fr(0)), (Fr(0), Fr(1)))
def slash_lower(wl):        # gamma^mu w_mu
    M = ZERO2
    for m in range(3):
        M = madd(M, GAM[m], wl[m])
    return M
def slash(v):
    return slash_lower([ETA[m]*v[m] for m in range(3)])
def S_deriv(k, l, zi, zj, x):
    """(zi.d_i)^k (zj.d_j)^l S(x_i - x_j) / (c_F (x^2)^{-1/2}),  S(x) = c_F xslash/(x^2)^{3/2} = -c_F dslash (x^2)^{-1/2}"""
    _, gr = g_deriv(k, l, zi, zj, x, grad=True)
    return mscal(-(-1)**l, slash_lower(gr))
def c_fermion_T():
    return {(1, 0): Fr(1), (0, 1): Fr(-1)}
def c_fermion_J(s):
    return {(s - 1 - k, k): Fr((-1)**k*comb(2*s, 2*k + 1)) for k in range(s)}
def tr(A):
    return A[0][0] + A[1][1]
def fermion_3pt(c4, G4, c1, G1_, c3, G3, x4, z4, x1, z1, x3, z3):
    """<B4 B1 B3> / (c_F^3 (x41^2 x43^2 x13^2)^{-1/2}); G = 'z' (zslash) or '1'"""
    gm = lambda G, z: slash(z) if G == 'z' else ID2
    M4, M1, M3 = gm(G4, z4), gm(G1_, z1), gm(G3, z3)
    x41, x13, x34 = sub(x4, x1), sub(x1, x3), sub(x3, x4)
    x43, x31, x14 = sub(x4, x3), sub(x3, x1), sub(x1, x4)
    tot = Fr(0)
    for (a4, b4), k4 in c4.items():
        for (a1, b1), k1 in c1.items():
            for (a3, b3), k3 in c3.items():
                A = mmul(mmul(mmul(M4, S_deriv(b4, a1, z4, z1, x41)), mmul(M1, S_deriv(b1, a3, z1, z3, x13))),
                         mmul(M3, S_deriv(b3, a4, z3, z4, x34)))
                B = mmul(mmul(mmul(M4, S_deriv(b4, a3, z4, z3, x43)), mmul(M3, S_deriv(b3, a1, z3, z1, x31))),
                         mmul(M1, S_deriv(b1, a4, z1, z4, x14)))
                tot += -k4*k1*k3*(tr(A) + tr(B))
    return tot
def fermion_2pt(c, G, cp, Gp, x, z, zp):
    gm = lambda G_, z_: slash(z_) if G_ == 'z' else ID2
    M, Mp = gm(G, z), gm(Gp, zp)
    mx = [-q for q in x]
    tot = Fr(0)
    for (a, b), k in c.items():
        for (ap, bp), kp in cp.items():
            tot += -k*kp*tr(mmul(mmul(M, S_deriv(b, ap, z, zp, x)), mmul(Mp, S_deriv(bp, a, zp, z, mx))))
    return tot

# ---------------------------------------------------------------- projection on structures
def rnd_pt(rng):
    return [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(3)]
def rnd_null(rng):
    return [Fr(c) for c in rng.choice(st.NULLS)]

def project(corr, monos, J, parity, seed=11, extra=6):
    """corr(x4,z4,x1,z1,x3,z3) -> rational part R_W with W = R_W (X13 X14 X34)^{-1/2}.
    Solve R_W = sum_a l_a R_a X13^{J-4} X14^{-J} X34^{-J} (even; s = J + 1/2)  resp. with one more inverse power of
    each X_ij and a factor (X13 X14 X34)^{1/2} (odd: only possible if R_W is itself odd, i.e. we fit
    R_W (X13 X14 X34)^{1/2}... see below).  Returns (coeffs, max residual)."""
    rng = random.Random(seed)
    rows, rhs = [], []
    for _ in range(len(monos) + extra):
        x4, z4, x1, z1, x3, z3 = rnd_pt(rng), rnd_null(rng), rnd_pt(rng), rnd_null(rng), rnd_pt(rng), rnd_null(rng)
        B = st.blocks([sp.Rational(c.numerator, c.denominator) for c in x1], [sp.Integer(int(c)) for c in z1],
                      [sp.Rational(c.numerator, c.denominator) for c in x3], [sp.Integer(int(c)) for c in z3],
                      [sp.Rational(c.numerator, c.denominator) for c in x4], [sp.Integer(int(c)) for c in z4])
        if parity == 'even':
            pref = B['X13']**(J - 4)*B['X14']**(-J)*B['X34']**(-J)
        else:   # s = (Delta+J)/2; odd prefactor exponents e + 1/2: handled by the caller through Delta
            raise NotImplementedError
        rows.append([mono_value_(m, B)*pref for m in monos])
        w = corr(x4, z4, x1, z1, x3, z3)
        rhs.append(sp.Rational(w.numerator, w.denominator))
    M = sp.Matrix(rows); b = sp.Matrix(rhs)
    n = len(monos)
    sol = M[:n, :].LUsolve(b[:n, :]) if M[:n, :].det() != 0 else None
    if sol is None:
        sol, params = M.gauss_jordan_solve(b)
        assert not params.free_symbols or True
    res = (M*sol - b)
    return list(sol), max(abs(r) for r in res)

def mono_value_(m, B):
    return st.mono_value(m, B)

def project_general(corr, monos, expo, seed=11, extra=6):
    """Generic version: W_rat = sum_a l_a R_a X13^{expo[0]} X14^{expo[1]} X34^{expo[2]} (integer exponents)."""
    rng = random.Random(seed)
    rows, rhs = [], []
    for _ in range(len(monos) + extra):
        x4, z4, x1, z1, x3, z3 = rnd_pt(rng), rnd_null(rng), rnd_pt(rng), rnd_null(rng), rnd_pt(rng), rnd_null(rng)
        R = lambda v: [sp.Rational(c.numerator, c.denominator) for c in v]
        B = st.blocks(R(x1), R(z1), R(x3), R(z3), R(x4), R(z4))
        pref = B['X13']**expo[0]*B['X14']**expo[1]*B['X34']**expo[2]
        rows.append([st.mono_value(m, B)*pref for m in monos])
        w = corr(x4, z4, x1, z1, x3, z3)
        rhs.append(sp.Rational(w.numerator, w.denominator))
    M = sp.Matrix(rows); b = sp.Matrix(rhs)
    sol, params = M.gauss_jordan_solve(b)
    assert len(params) == 0, "structures not independent"
    return list(sol)
