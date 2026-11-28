"""First-order test of the crossing inequality (Lemma X) and of Lemma N3 in a
two-sided BTZ black hole (r_h = 1, pi_1 = Z) perturbed by a thin null shell
thrown in from the right boundary at boundary time v0 (Vaidya-BTZ).

Exact Vaidya metric in ingoing EF coordinates of the right side:
    ds^2 = -(r^2 - M(v)) dv^2 + 2 dv dr + r^2 dphi^2,  M = 1 + dM * theta(v - v0).
So g = g_BTZ + dM * theta(v-v0) dv^2 EXACTLY, and the renormalized length of the
extremal geodesic in a fixed homotopy class is
    D = D0 + dM * D1 + O(dM^2),   D1 = 1/2 int theta(v-v0) (dv/ds)^2 ds
along the vacuum BTZ geodesic.  NEC  <=>  dM >= 0.

Embedding: -T1^2 - T2^2 + X1^2 + X2^2 = -1.  Right boundary point (t, phi):
P = (T1,T2,X1,X2) = (cosh phi, sinh t, cosh t, sinh phi); left: X1 -> -X1.
Ingoing EF time (valid where V = X1 + T2 > 0): v = log(V/(1+r)), r^2 = T1^2 - X2^2.
Universal cover: phi in R (deck group phi -> phi + 2 pi).
Renormalized vacuum length: D0 = log(2 * (-P.Q)).

usage: python3 -I vaidya_btz_first_order.py [seed]
"""
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

V0 = 0.3          # shell injection time on the right boundary


def bpt(side, t, phi):
    s = 1.0 if side == 'R' else -1.0
    return np.array([np.cosh(phi), np.sinh(t), s * np.cosh(t), np.sinh(phi)])


def dot(a, b):
    return -a[0] * b[0] - a[1] * b[1] + a[2] * b[2] + a[3] * b[3]


def D0(P, Q):
    return np.log(2.0 * (-dot(P, Q)))


def D1(P, Q, v0=V0):
    n = np.sqrt(-2.0 * dot(P, Q))

    def x(l):
        return (np.exp(l) * P + np.exp(-l) * Q) / n

    def xd(l):
        return (np.exp(l) * P - np.exp(-l) * Q) / n

    def vfun(l):
        X = x(l)
        V = X[2] + X[1]
        r2 = X[0] ** 2 - X[3] ** 2
        if V <= 0 or r2 <= 0:
            return -np.inf
        return np.log(V / (1.0 + np.sqrt(r2)))

    def vdot2(l):
        X, Xd = x(l), xd(l)
        V, Vd = X[2] + X[1], Xd[2] + Xd[1]
        r = np.sqrt(X[0] ** 2 - X[3] ** 2)
        rd = (X[0] * Xd[0] - X[3] * Xd[3]) / r
        return (Vd / V - rd / (1.0 + r)) ** 2

    # locate the set {v > v0} on a grid, refine its endpoints
    L = 40.0
    grid = np.linspace(-L, L, 4001)
    vals = np.array([vfun(l) for l in grid]) - v0
    tot = 0.0
    inside = vals > 0
    i = 0
    N = len(grid)
    while i < N:
        if inside[i]:
            j = i
            while j + 1 < N and inside[j + 1]:
                j += 1
            a = grid[i] if i == 0 else brentq(lambda l: vfun(l) - v0, grid[i - 1], grid[i])
            b = grid[j] if j == N - 1 else brentq(lambda l: vfun(l) - v0, grid[j], grid[j + 1])
            tot += quad(vdot2, a, b, limit=200, epsabs=1e-13, epsrel=1e-11)[0]
            i = j + 1
        else:
            i += 1
    return 0.5 * tot


def slice_times(side, phi, pars):
    # boundary Cauchy slice: periodic, |dt/dphi| < 1
    if side == 'R':
        a0, a1, c1 = pars['R']
        return a0 + a1 * np.sin(phi + c1)
    a0, a1, c1 = pars['L']
    return a0 + a1 * np.sin(2 * phi + c1)


def key(pt):
    side, phi = pt
    return (0, phi) if side == 'R' else (1, -phi)


def point(pt, pars):
    side, phi = pt
    return bpt(side, slice_times(side, phi, pars), phi)


def rect(pts, pars):
    """pts cyclically ordered x1<x2<y1<y2; returns (B0, B1) with
    nu = 1/2 [D(x1,y1)+D(x2,y2)-D(x1,y2)-D(x2,y1)] = B0 + dM*B1."""
    x1, x2, y1, y2 = [point(p, pars) for p in pts]
    # B0 via a single log of the cross-ratio (no cancellation)
    B0 = 0.5 * np.log((dot(x1, y1) * dot(x2, y2)) / (dot(x1, y2) * dot(x2, y1)))
    B1 = 0.5 * (D1(x1, y1) + D1(x2, y2) - D1(x1, y2) - D1(x2, y1))
    return B0, B1


def random_quad(rng, scale):
    pts = []
    for _ in range(4):
        side = 'R' if rng.random() < 0.5 else 'L'
        pts.append((side, rng.normal(0, scale)))
    pts.sort(key=key)
    return pts


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    rng = np.random.default_rng(seed)
    pars = {'R': (0.6, 0.4, 0.3), 'L': (0.4, 0.3, 1.1)}
    # ---- T1: random quadruples
    res = []
    for k in range(400):
        pts = random_quad(rng, 2.5)
        if len(set((p[0], round(p[1], 12)) for p in pts)) < 4:
            continue
        B0, B1 = rect(pts, pars)
        res.append((B0, B1))
    res = np.array(res)
    print('T1 random: n=%d  min B0=%.3e  (B0>0 all: %s)' % (len(res), res[:, 0].min(), (res[:, 0] > 0).all()))
    print('   max(-B1/B0) [NEC sign threat] = %.4f ; max(B1/B0) [anti-NEC threat] = %.4f'
          % ((-res[:, 1] / res[:, 0]).max(), (res[:, 1] / res[:, 0]).max()))
    # ---- T1b: shrinking rectangles straddling a tangency locus (R-L pairs)
    # scan D1 along a right-left family to find where the geodesic starts touching the shell
    pR = 0.0
    phis = np.linspace(-3, 3, 241)
    d1 = np.array([D1(point(('R', pR), pars), point(('L', p), pars)) for p in phis])
    on = np.where(d1 > 1e-12)[0]
    print('T1b: R-L family, D1>0 for phi_L in [%.3f, %.3f] (onset = tangency)' %
          (phis[on[0]], phis[on[-1]]) if len(on) else 'T1b: no crossing in family')
    if len(on) and on[0] > 0:
        # refine onset
        f = lambda p: D1(point(('R', pR), pars), point(('L', p), pars)) - 1e-14
        lo, hi = phis[on[0] - 1], phis[on[0]]
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        p_star = hi
        print('     onset phi_L* = %.10f' % p_star)
        for h in [0.2, 0.05, 0.0125, 0.003, 0.0008]:
            worst_nec, worst_anti = 0, 0
            for off in np.linspace(-1.5 * h, 1.5 * h, 7):
                # x1<x2 on right near pR, y1<y2 on left near p_star (left order: phi decreasing)
                pts = [('R', pR - h / 2), ('R', pR + h / 2), ('L', p_star + off + h / 2), ('L', p_star + off - h / 2)]
                B0, B1 = rect(pts, pars)
                worst_nec = max(worst_nec, -B1 / B0)
                worst_anti = max(worst_anti, B1 / B0)
            print('     h=%.4f  max(-B1/B0)=%9.4f   max(B1/B0)=%9.4f' % (h, worst_nec, worst_anti))
    # ---- T2: Lemma N3, d_n -> L(h) = 2 pi, first-order shift -> 0
    x = ('R', 0.2)
    y = ('L', -0.4)
    Dn = []
    for n in range(0, 9):
        P = point(x, pars)
        Q = bpt('L', slice_times('L', y[1], pars), y[1] + 2 * np.pi * n)  # h^n y (deck shift)
        Dn.append((D0(P, Q), D1(P, Q)))
    Dn = np.array(Dn)
    d0 = np.diff(Dn[:, 0])
    d1n = np.diff(Dn[:, 1])
    print('T2: d_n^(0) - 2pi =', ' '.join('%.2e' % (a - 2 * np.pi) for a in d0))
    print('    d_n^(1)       =', ' '.join('%.2e' % a for a in d1n))
    print('    d_n^(0) nondecreasing:', bool(np.all(np.diff(d0) >= -1e-12)))


if __name__ == '__main__':
    main()
