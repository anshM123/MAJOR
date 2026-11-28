"""Fast double-precision flat thermal bubble Pi0_n(u) and u-derivatives (beta = 1, m0 = 2 log phi).
Vectorised panel Gauss-Legendre in q; validated against flat.py (mpmath) in check_fastflat.py.
"""
import numpy as np

M0 = 2 * np.log((1 + np.sqrt(5)) / 2)
_x, _w = np.polynomial.legendre.leggauss(48)


def _panels(edges):
    Q, W = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        Q.append(0.5 * (b - a) * _x + 0.5 * (b + a))
        W.append(0.5 * (b - a) * _w)
    return np.concatenate(Q), np.concatenate(W)


def pi_vac(u, w, m=M0):
    """returns (Pi, dPi/du, d2Pi/du2) of atan(p/2m)/(4 pi p), p^2 = u + w^2 (p > 0 assumed or series)."""
    u = np.asarray(u, float)
    p2 = u + w * w
    p = np.sqrt(p2)
    x = p / (2 * m)
    out = []
    small = x < 0.3
    # series for small x in X = x^2:  atan(x)/x = sum (-1)^k X^k/(2k+1)
    X = x * x
    ks = np.arange(0, 40)
    c = (-1.0) ** ks / (2 * ks + 1)
    s0 = sum(c[k] * X ** k for k in ks)
    s1 = sum(c[k] * k * X ** (k - 1) for k in ks[1:])
    s2 = sum(c[k] * k * (k - 1) * X ** (k - 2) for k in ks[2:])
    A = np.arctan(x)
    with np.errstate(divide='ignore', invalid='ignore'):
        D = 1 / (1 + x * x) / (2 * m)
        DD = -2 * x / (1 + x * x) ** 2 / (2 * m) ** 2
        g0 = A / (4 * np.pi * p)
        g1 = (D / p - A / p2) / (4 * np.pi)
        g2 = (DD / p - 2 * D / p2 + 2 * A / (p2 * p)) / (4 * np.pi)
        f1 = g1 / (2 * p)
        f2 = g2 / (4 * p2) - g1 / (4 * p2 * p)
    pre = 1 / (8 * np.pi * m)
    f0 = np.where(small, pre * s0, g0)
    f1 = np.where(small, pre * s1 / (4 * m * m), f1)
    f2 = np.where(small, pre * s2 / (4 * m * m) ** 2, f2)
    return f0, f1, f2


def pi_th(n, u, m=M0):
    """thermal part and u-derivatives at Matsubara n, array u."""
    u = np.atleast_1d(np.asarray(u, float))
    if n == 0:
        res = np.empty((3, len(u)))
        for i, uu in enumerate(u):
            sc = 2 / np.sqrt(uu + 1.0)
            ed = np.concatenate([[0], sc * 2.0 ** np.arange(-1, 9)])
            ed = np.unique(np.concatenate([ed[ed < np.pi / 2], np.pi / 2 * np.arange(1, 5) / 4]))
            ed = ed[ed <= np.pi / 2]
            th, wt = _panels(ed)
            s2 = np.sin(th) ** 2
            E = np.sqrt(uu * s2 / 4 + m * m)
            with np.errstate(over='ignore', invalid='ignore'):
                qe = np.exp(-E)
                om = -np.expm1(-E)
                g = qe / om / E
                r1 = qe / om ** 2              # = -n'(E)
                r2 = qe * (1 + qe) / om ** 3   # = n''(E)
                gE = -r1 / E - g / E
                gEE = r2 / E + 2 * r1 / E ** 2 + 2 * g / E ** 2
            big = E > 600
            g = np.where(big, 0, g); gE = np.where(big, 0, gE); gEE = np.where(big, 0, gEE)
            Eu = s2 / (8 * E)
            Euu = -s2 / (8 * E * E) * Eu
            sw = np.sin(th) * wt / (4 * np.pi)
            res[:, i] = (g @ sw, (gE * Eu) @ sw, (gEE * Eu * Eu + gE * Euu) @ sw)
        return res[0], res[1], res[2]
    w = 2 * np.pi * n
    edges = np.array([0, 0.5, 1, 2, 3, 4.5, 6, 8, 11, 15, 20, 27, 36, 48, 64])
    # add panel edges near the quasi-singular point q0 = sqrt(u)/2 for each u (done per u-chunk below)
    out0 = np.empty_like(u); out1 = np.empty_like(u); out2 = np.empty_like(u)
    for i, uu in enumerate(u):
        q0 = np.sqrt(uu) / 2
        ed = edges
        if 0.5 < q0 < 64:
            extra = q0 + np.array([-2, -1, -0.5, -0.2, 0, 0.2, 0.5, 1, 2]) * max(1.0, w * 0.5)
            ed = np.unique(np.concatenate([edges, extra[(extra > 0) & (extra < 64)]]))
        q, wt = _panels(ed)
        E = np.sqrt(q * q + m * m)
        a = -(uu + w * w) + 2j * w * E
        b2 = 4 * uu * q * q
        S = a * np.sqrt(1 - b2 / (a * a))
        wq = q / (np.expm1(E) * E) * wt
        out0[i] = -(np.real(1 / S) @ wq) / np.pi
        out1[i] = -(np.real((a + 2 * q * q) / S ** 3) @ wq) / np.pi
        out2[i] = -(np.real(-1 / S ** 3 + 3 * (a + 2 * q * q) ** 2 / S ** 5) @ wq) / np.pi
    return out0, out1, out2


def Pi0(n, u):
    w = 2 * np.pi * n
    v = pi_vac(u, w)
    t = pi_th(n, u)
    return v[0] + t[0], v[1] + t[1], v[2] + t[2]
