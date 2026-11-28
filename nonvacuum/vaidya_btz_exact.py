"""EXACT thin-shell Vaidya-BTZ test of the crossing inequality (Lemma X) on the
universal cover, and of Lemma N3 (d_n -> L(h)), for NEC (dM>0) and NEC-violating
(dM<0) shells.

Two-sided BTZ of mass M1 = 1; a null shell enters from the RIGHT boundary at
boundary time v0, after which the mass is M2 = M1 + dM (ingoing EF chart:
ds^2 = -(r^2-M(v))dv^2 + 2dvdr + r^2 dphi^2).  The left exterior is pure BTZ(M1).

A right-left geodesic from P (right boundary, time tR > v0) to Q (left boundary)
is a BTZ(M2) geodesic from P to a break point c on the shell {v=v0}, then a BTZ(M1)
geodesic from c to Q.  The true geodesic is the critical point of the total
renormalized length over the break point c = (r_c, phi_c) (Fermat; equivalent to
the refraction law  dot v, dot phi continuous, dot r jumps by (M_before-M_after) dot v/2).
Each half has renormalized length log(2(-X.B)) in the embedding of its own BTZ,
B the boundary null vector normalised by X/r -> B.

All four lengths in a rectangle  x1<x2 (right) , y1<y2 (left)  are right-left,
so nu(rectangle) = 1/2[D(x1,y1)+D(x2,y2)-D(x1,y2)-D(x2,y1)] is computed exactly.

usage: python3 -I vaidya_btz_exact.py
"""
import numpy as np
from scipy.optimize import root

V0 = 0.0


def emb(M, v, r, phi):
    m = np.sqrt(M)
    T1 = (r / m) * np.cosh(m * phi)
    X2 = (r / m) * np.sinh(m * phi)
    V = (r + m) * np.exp(m * v) / m
    U = (r * r / (m * m) - 1.0) / V
    return np.array([T1, 0.5 * (V - U), 0.5 * (V + U), X2])  # (T1,T2,X1,X2)


def bnull(M, side, t, phi):
    m = np.sqrt(M)
    s = 1.0 if side == 'R' else -1.0
    return np.array([np.cosh(m * phi), np.sinh(m * t), s * np.cosh(m * t), np.sinh(m * phi)]) / m


def dot(a, b):
    return -a[0] * b[0] - a[1] * b[1] + a[2] * b[2] + a[3] * b[3]


def vcoord(M, X):
    m = np.sqrt(M)
    V = X[2] + X[1]
    r2 = X[0] ** 2 - X[3] ** 2
    if V <= 0 or r2 <= 0:
        return -np.inf
    r = m * np.sqrt(r2)
    return np.log(V * m / (r + m)) / m


def D0(tR, pR, tL, pL):
    P, Q = bnull(1.0, 'R', tR, pR), bnull(1.0, 'L', tL, pL)
    return np.log(2.0 * (-dot(P, Q)))


def vac_cross(tR, pR, tL, pL, v0=V0):
    """crossing point of the vacuum geodesic with {v = v0} (initial guess)."""
    P, Q = bnull(1.0, 'R', tR, pR), bnull(1.0, 'L', tL, pL)
    n = np.sqrt(-2 * dot(P, Q))
    x = lambda l: (np.exp(l) * P + np.exp(-l) * Q) / n
    ls = np.linspace(40, -40, 8001)
    vs = np.array([vcoord(1.0, x(l)) for l in ls])
    k = np.where(vs < v0)[0][0]
    l0 = ls[k]
    X = x(l0)
    r = np.sqrt(X[0] ** 2 - X[3] ** 2)
    phi = np.arctanh(X[3] / X[0])
    return r, phi


def F(c, M2, P2, Q1, v0=V0):
    r, phi = c
    X2 = emb(M2, v0, r, phi)
    X1 = emb(1.0, v0, r, phi)
    a, b = -dot(X2, P2), -dot(X1, Q1)
    if a <= 0 or b <= 0:
        return np.nan
    return np.log(2 * a) + np.log(2 * b)


def grad(c, *args, h=1e-6):
    g = np.zeros(2)
    for i in range(2):
        e = np.zeros(2); e[i] = h
        g[i] = (F(c + e, *args) - F(c - e, *args)) / (2 * h)
    return g


def D_exact(dM, tR, pR, tL, pL, guess=None, v0=V0):
    M2 = 1.0 + dM
    P2 = bnull(M2, 'R', tR, pR)
    Q1 = bnull(1.0, 'L', tL, pL)
    c0 = np.array(guess if guess is not None else vac_cross(tR, pR, tL, pL, v0))
    sol = root(grad, c0, args=(M2, P2, Q1, v0), method='hybr', tol=1e-13)
    c = sol.x
    val = F(c, M2, P2, Q1, v0)
    # validity: M2 part stays in v>=v0, M1 part in v<=v0
    X2 = emb(M2, v0, *c); X1 = emb(1.0, v0, *c)
    ok = bool(np.isfinite(val) and np.linalg.norm(grad(c, M2, P2, Q1, v0)) < 1e-7)
    for s in np.linspace(0.05, 12, 40):
        y2 = np.exp(-s) * X2 + np.sinh(s) * P2 / (-dot(X2, P2))
        y1 = np.exp(-s) * X1 + np.sinh(s) * Q1 / (-dot(X1, Q1))
        if vcoord(M2, y2) < v0 - 1e-9:
            ok = False
        v1 = vcoord(1.0, y1)
        if np.isfinite(v1) and v1 > v0 + 1e-9:
            ok = False
    return val, c, ok


def check_vacuum_limit():
    tR, pR, tL, pL = 0.5, 0.1, 0.3, -0.7
    d, c, ok = D_exact(0.0, tR, pR, tL, pL)
    print('vacuum limit: D_exact(dM=0) - D0 = %.2e  (ok=%s)' % (d - D0(tR, pR, tL, pL), ok))
    eps = 1e-5
    dp, _, _ = D_exact(eps, tR, pR, tL, pL)
    print('  dD/ddM at 0 (finite diff) = %.6f' % ((dp - d) / eps))


def rect_scan(dM, tR, tL, rng, n=300, scale=2.0):
    worst = np.inf
    nbad = 0
    nok = 0
    for _ in range(n):
        a, b = np.sort(rng.normal(0, scale, 2))
        y = np.sort(rng.normal(0, scale, 2))[::-1]  # left points: cyclic order = phi decreasing
        pts = [(a, y[0]), (b, y[1]), (a, y[1]), (b, y[0])]
        vals = []
        good = True
        for (pr, pl) in pts:
            d, c, ok = D_exact(dM, tR, pr, tL, pl)
            good &= ok
            vals.append(d)
        if not good:
            continue
        nok += 1
        nu = 0.5 * (vals[0] + vals[1] - vals[2] - vals[3])
        nu0 = 0.5 * (D0(tR, a, tL, y[0]) + D0(tR, b, tL, y[1]) - D0(tR, a, tL, y[1]) - D0(tR, b, tL, y[0]))
        worst = min(worst, nu / nu0)
        nbad += nu < -1e-11
    return nok, nbad, worst


def convexity_scan(dM, tR, tL, pR=0.0, span=4.0, N=161):
    """for constant-time slices nu>=0 on R-L rectangles <=> f(phi_L) = D(pR, phi_L) convex."""
    ps = np.linspace(-span, span, N)
    f = []
    guess = None
    for p in ps:
        d, c, ok = D_exact(dM, tR, pR, tL, p, guess)
        f.append(d if ok else np.nan)
        guess = c if ok else None
    f = np.array(f)
    h = ps[1] - ps[0]
    f2 = (f[2:] - 2 * f[1:-1] + f[:-2]) / h ** 2
    return np.nanmin(f2), np.nanmax(f2), np.isnan(f).sum()


def N3_test(dM, tR, tL, pR=0.0, pL=0.4, nmax=6):
    Ds = []
    guess = None
    for n in range(nmax + 1):
        d, c, ok = D_exact(dM, tR, pR, tL, pL + 2 * np.pi * n, guess)
        Ds.append(d if ok else np.nan)
        guess = None
    return np.diff(Ds)


def main():
    check_vacuum_limit()
    rng = np.random.default_rng(1)
    tR, tL = 0.6, 0.2   # right endpoints after the shell (v0 = 0), left before
    for dM in [2.0, 0.5, 0.0, -0.3, -0.6, -0.85]:
        mn, mx, nn = convexity_scan(dM, tR, tL)
        nok, nbad, worst = rect_scan(dM, tR, tL, rng, n=120)
        d = N3_test(dM, tR, tL)
        print('dM=%+.2f  min f\'\'=%+.4f max f\'\'=%+.4f (nan %d) | random rects ok=%d  nu<0: %d  min nu/nu0=%.4f'
              % (dM, mn, mx, nn, nok, nbad, worst))
        print('          d_n - 2pi =', ' '.join('%+.2e' % (x - 2 * np.pi) for x in d))


if __name__ == '__main__':
    main()
