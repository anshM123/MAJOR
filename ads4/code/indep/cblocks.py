"""Continuous-spin d = 3 conformal blocks with external shift (a, b) = (0, b), Delta = 1 (or general), spin l = -1/2 + i nu,
in Lorentzian radial variables rho = r e^beta, rhobar = r e^-beta, eta = cosh beta, z = 4 rho/(1+rho)^2.
Own derivation (see INDEP_CHECK.md, Target A methods):

  Dolan-Osborn 1/2-Casimir, a = 0:
     z^2(1-z)d_z^2 + zb^2(1-zb)d_zb^2 - (b+1)(z^2 d_z + zb^2 d_zb) + (d-2) z zb/(z-zb) ((1-z)d_z - (1-zb)d_zb),
     eigenvalue C/2, C = Delta(Delta-d) + l(l+d-2).
  In (r, eta), d = 3, multiplied by 2 (1-r^2) Den, Den = (1+r^2)^2 - 4 r^2 eta^2:
     (1-r^2) Den (D_r^2 + Lam - C) G - [2(1-r^2)^2((1+r^2) + 4 b r eta) + Den (1+r^2)] D_r G
        - 8 r (1-r^2) [r eta + b(1+r^2)] (eta^2-1) d_eta G = 0,
     D_r = r d_r, Lam = d_eta (eta^2-1) d_eta  (Legendre: Lam P_mu = mu(mu+1) P_mu).
  G = r^Delta sum_N r^N sum_m B[N,m] P_{l+m}(eta), B[0,0] = 1  (Ghat := G / r^Delta -> P_l as r -> 0).
"""
import numpy as np
import mpmath as mp

# ---------------------------------------------------------------- operator polynomials
def _poly_terms(b):
    """Return lists of (k, j, coeff) for C1 (mult. D^2-C+Lam), C2 (mult. D), C3 (mult. D1=(eta^2-1)d_eta):
    coefficient sum_k,j coeff r^k eta^j."""
    import sympy as sp
    r, e, bb = sp.symbols('r e bb')
    Den = (1 + r**2)**2 - 4 * r**2 * e**2
    C1 = sp.expand((1 - r**2) * Den)
    C2 = sp.expand(-(2 * (1 - r**2)**2 * ((1 + r**2) + 4 * bb * r * e) + Den * (1 + r**2)))
    C3 = sp.expand(-8 * r * (1 - r**2) * (r * e + bb * (1 + r**2)))
    out = []
    for Cx in (C1, C2, C3):
        P = sp.Poly(Cx, r, e)
        out.append([(int(k), int(j), complex(sp.N(c.subs(bb, b)))) for (k, j), c in zip(P.monoms(), P.coeffs())])
    return out

def recursion(nu, b, Nmax, Delta=1):
    """B[N][m] (dict m -> complex) for N = 0..Nmax."""
    l = -0.5 + 1j * nu
    C = Delta * (Delta - 3) + l * (l + 1)
    T1, T2, T3 = _poly_terms(b)
    B = [dict() for _ in range(Nmax + 1)]
    B[0][0] = 1.0 + 0j

    def eta_mult(g):
        out = {}
        for m, c in g.items():
            mu = l + m
            out[m + 1] = out.get(m + 1, 0) + c * (mu + 1) / (2 * mu + 1)
            out[m - 1] = out.get(m - 1, 0) + c * mu / (2 * mu + 1)
        return out

    def D1(g):
        out = {}
        for m, c in g.items():
            mu = l + m
            f = mu * (mu + 1) / (2 * mu + 1)
            out[m + 1] = out.get(m + 1, 0) + c * f
            out[m - 1] = out.get(m - 1, 0) - c * f
        return out

    def add(acc, g, fac):
        for m, c in g.items():
            acc[m] = acc.get(m, 0) + fac * c

    for N in range(1, Nmax + 1):
        src = {}
        for (terms, kind) in ((T1, 1), (T2, 2), (T3, 3)):
            for (k, j, cf) in terms:
                if k == 0 or k > N:
                    continue
                g = B[N - k]
                if not g:
                    continue
                D = Delta + N - k
                if kind == 1:
                    h = {m: c * (D * D - C + (l + m) * (l + m + 1)) for m, c in g.items()}
                elif kind == 2:
                    h = {m: c * D for m, c in g.items()}
                else:
                    h = D1(g)
                for _ in range(j):
                    h = eta_mult(h)
                add(src, h, cf)
        D = Delta + N
        for m, c in src.items():
            E = D * D - 3 * D + (l + m) * (l + m + 1) - C   # Op_0 eigenvalue (C1, C2 at r^0: 1 and -3)
            if abs(c) == 0:
                continue
            if abs(E) < 1e-300:
                raise ZeroDivisionError((N, m))
            B[N][m] = -c / E
    return B

# ---------------------------------------------------------------- Legendre functions P_{l+m}(cosh beta)
def legendre_table(nu, eta, Mmax):
    """P_{l+m}(eta) for m = -Mmax..Mmax, l = -1/2 + i nu, eta > 1 (upward recursion; P_{l-m} = conj P_{l+m})."""
    l = mp.mpc(-0.5, nu)
    p0 = complex(mp.legenp(l, 0, eta, type=3))
    p1 = complex(mp.legenp(l + 1, 0, eta, type=3))
    P = {0: p0, 1: p1}
    lc = complex(l)
    for m in range(1, Mmax):
        mu = lc + m
        P[m + 1] = ((2 * mu + 1) * eta * P[m] - mu * P[m - 1]) / (mu + 1)
    for m in range(1, Mmax + 1):
        P[-m] = np.conj(P[m])
    return P

def Ghat(nu, b, r, eta, Nmax=200, B=None, P=None):
    if B is None:
        B = recursion(nu, b, Nmax)
    if P is None:
        P = legendre_table(nu, eta, Nmax + 1)
    s = 0j
    rn = 1.0
    for N in range(len(B)):
        t = 0j
        for m, c in B[N].items():
            t += c * P[m]
        s += rn * t
        rn *= r
    return s

# ---------------------------------------------------------------- checks of the operator
def check_operators(log=print):
    mp.mp.dps = 30
    # (1) DO operator in (z, zb) reproduces the d = 4 closed-form block with (a, b) = (0, b)
    bv = mp.mpf('0.7'); a = 0
    Dl, ll = mp.mpf('2.3'), 2
    def k(beta, x): return x**(beta / 2) * mp.hyp2f1(beta / 2 + a, beta / 2 + bv, beta, x)
    def G4(z, zb): return z * zb / (z - zb) * (k(Dl + ll, z) * k(Dl - ll - 2, zb) - k(Dl + ll, zb) * k(Dl - ll - 2, z))
    def DO(f, z, zb, d):
        Gz = mp.diff(lambda t: f(t, zb), z); Gzz = mp.diff(lambda t: f(t, zb), z, 2)
        Gw = mp.diff(lambda t: f(z, t), zb); Gww = mp.diff(lambda t: f(z, t), zb, 2)
        return (z**2 * (1 - z) * Gzz + zb**2 * (1 - zb) * Gww - (a + bv + 1) * (z**2 * Gz + zb**2 * Gw)
                - a * bv * (z + zb) * f(z, zb) + (d - 2) * z * zb / (z - zb) * ((1 - z) * Gz - (1 - zb) * Gw))
    z, zb = mp.mpf('0.31'), mp.mpf('0.12')
    lhs = DO(G4, z, zb, 4); rhs = (Dl * (Dl - 4) + ll * (ll + 2)) / 2 * G4(z, zb)
    log(f'  DO operator vs d=4 closed-form block (b=0.7): {mp.nstr(lhs, 18)} vs {mp.nstr(rhs, 18)}  rel {mp.nstr(abs(lhs/rhs-1), 3)}')
    # (2) radial form vs DO operator on a test function, d = 3
    def f_zz(z, zb): return mp.exp(z) * mp.cos(zb) + z**2 * zb**3 + mp.log(1 + z * zb)
    def f_re(r, e):
        s = mp.sqrt(e * e - 1); rho, rhb = r * (e + s), r * (e - s)
        return f_zz(4 * rho / (1 + rho)**2, 4 * rhb / (1 + rhb)**2)
    C = mp.mpf(0)   # test with C = 0: compare 2(1-r^2)Den * DO[f] with the radial expression (C-term absent)
    for (r, e) in [(mp.mpf('0.2'), mp.mpf('1.3')), (mp.mpf('0.35'), mp.mpf('2.1'))]:
        Den = (1 + r**2)**2 - 4 * r**2 * e**2
        Dr = r * mp.diff(lambda t: f_re(t, e), r)
        Dr2 = r * mp.diff(lambda t: t * mp.diff(lambda u: f_re(u, e), t), r)
        fe = mp.diff(lambda t: f_re(r, t), e); fee = mp.diff(lambda t: f_re(r, t), e, 2)
        Lam = (e * e - 1) * fee + 2 * e * fe
        D1 = (e * e - 1) * fe
        rad = ((1 - r**2) * Den * (Dr2 + Lam) - (2 * (1 - r**2)**2 * ((1 + r**2) + 4 * bv * r * e) + Den * (1 + r**2)) * Dr
               - 8 * r * (1 - r**2) * (r * e + bv * (1 + r**2)) * D1)
        s = mp.sqrt(e * e - 1); rho, rhb = r * (e + s), r * (e - s)
        do = 2 * (1 - r**2) * Den * DO(f_zz, 4 * rho / (1 + rho)**2, 4 * rhb / (1 + rhb)**2, 3)
        log(f'  radial operator vs DO (d=3, b=0.7) at r={r}, eta={e}: {mp.nstr(rad, 18)} vs {mp.nstr(do, 18)}  rel {mp.nstr(abs(rad/do-1), 3)}')

if __name__ == '__main__':
    check_operators()
    import time
    t = time.time(); B = recursion(1.3, 2, 120); print('recursion time', time.time() - t)
    print('B[1]', B[1]); print('B[2]', B[2])

# ---------------------------------------------------------------- vectorised over nu (same recursion)
def recursion_vec(nus, b, Nmax, Delta=1):
    nus = np.asarray(nus, float)
    l = -0.5 + 1j * nus
    C = Delta * (Delta - 3) + l * (l + 1)
    T1, T2, T3 = _poly_terms(b)
    B = [dict() for _ in range(Nmax + 1)]
    B[0][0] = np.ones_like(l)
    def eta_mult(g):
        out = {}
        for m, c in g.items():
            mu = l + m
            out[m + 1] = out.get(m + 1, 0) + c * (mu + 1) / (2 * mu + 1)
            out[m - 1] = out.get(m - 1, 0) + c * mu / (2 * mu + 1)
        return out
    def D1(g):
        out = {}
        for m, c in g.items():
            mu = l + m; f = mu * (mu + 1) / (2 * mu + 1)
            out[m + 1] = out.get(m + 1, 0) + c * f; out[m - 1] = out.get(m - 1, 0) - c * f
        return out
    for N in range(1, Nmax + 1):
        src = {}
        for (terms, kind) in ((T1, 1), (T2, 2), (T3, 3)):
            for (k, j, cf) in terms:
                if k == 0 or k > N or not B[N - k]:
                    continue
                g = B[N - k]; D = Delta + N - k
                if kind == 1:   h = {m: c * (D * D - C + (l + m) * (l + m + 1)) for m, c in g.items()}
                elif kind == 2: h = {m: c * D for m, c in g.items()}
                else:           h = D1(g)
                for _ in range(j): h = eta_mult(h)
                for m, c in h.items(): src[m] = src.get(m, 0) + cf * c
        D = Delta + N
        for m, c in src.items():
            if np.all(c == 0): continue
            E = D * D - 3 * D + (l + m) * (l + m + 1) - C
            B[N][m] = -c / E
    return B

def scaled_legendre_vec(nus, r, eta, Mmax):
    """Q_m = r^|m| P_{l+m}(eta), m = -Mmax..Mmax, vectorised over nu (avoids overflow at large eta)."""
    nus = np.asarray(nus, float)
    l = -0.5 + 1j * nus
    P0 = np.array([complex(mp.legenp(mp.mpc(-0.5, x), 0, eta, type=3)) for x in nus])
    P1 = np.array([complex(mp.legenp(mp.mpc(0.5, x), 0, eta, type=3)) for x in nus])
    Q = {0: P0, 1: r * P1}
    for m in range(1, Mmax):
        mu = l + m
        Q[m + 1] = ((2 * mu + 1) * (r * eta) * Q[m] - mu * r * r * Q[m - 1]) / (mu + 1)
    for m in range(1, Mmax + 1):
        Q[-m] = np.conj(Q[m])
    return Q

def Ghat_vec(B, Q, r):
    s = 0
    for N in range(len(B)):
        for m, c in B[N].items():
            s = s + c * Q[m] * r**(N - abs(m))
    return s
