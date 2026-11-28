"""n = 0 Matsubara sector:  S0(r) = sum_J (2J+1) [log P(u_J) - log Pi_vac(0,J)],  u_J = J(J+1)/r^2,
P(u) = Pi0_0(u) - (K/3) u Pi0_0''(u)   (curvature-expanded exact S^2 bubble, error O(K^2)),
Pi_vac(0,J) = (r/16) Gamma((J+1)/2)^2 / Gamma((J+2)/2)^2.
Prediction: S0(r) = A0 r^2 + (1/6) log r + B0 + o(1),  A0 = int_0^inf du log(8 sqrt(u) Pi0_0(u)).
Direct sum for J <= Jm, Euler-Maclaurin (midpoint form) for the tail with the exact continuous-J summand.
"""
import sys, numpy as np
from scipy.special import gammaln
import fastflat as ff

xg, wg = np.polynomial.legendre.leggauss(48)


from mpmath import bernpoly, mpf as _mpf
_YC = [(k, float(-4 * bernpoly(k + 1, _mpf(1) / 4) / (k * (k + 1)) * 2 ** k)) for k in range(2, 30, 2)]


def Y_of_l(l):
    """Y(l) = log(l V(l)), V = (1/2) Gamma((l+1/2)/2)^2/Gamma((l+3/2)/2)^2; series for l > 40."""
    l = np.asarray(l, float)
    big = l > 40
    lb = np.where(big, l, 50.0)
    ser = sum(c * lb ** (-k) for k, c in _YC)
    ls = np.where(big, 1.0, l)
    ex = np.log(ls) - np.log(2.0) + 2 * (gammaln((ls + 0.5) / 2) - gammaln((ls + 1.5) / 2))
    return np.where(big, ser, ex)


def F(J, r):
    """(2J+1)[log P(u_J) - log Pi_vac(0,J)] written as 2l[log(8 sqrt(u) P) + (1/2)log(l^2/(l^2-1/4)) - Y(l)]
    with every piece computed without cancellation; J = 0 handled exactly."""
    J = np.asarray(J, float)
    K = 1.0 / r ** 2
    l = J + 0.5
    u = J * (J + 1) * K
    out = np.empty_like(J)
    z = J == 0
    if z.any():
        P00 = ff.Pi0(0, np.array([0.0]))[0][0]
        out[z] = np.log(P00) - (np.log(r / 16.0) + 2 * (gammaln(0.5) - gammaln(1.0)))
    nz = ~z
    if nz.any():
        uu = u[nz]
        th, th1, th2 = ff.pi_th(0, uu)
        v0, v1, v2 = ff.pi_vac(uu, 0.0)
        sq = np.sqrt(uu)
        eps = -(2 / np.pi) * np.arctan(2 * ff.M0 / sq) + 8 * sq * th
        corr = -(K / 3) * uu * 8 * sq * (v2 + th2)
        X = np.log1p(eps + corr) - 0.5 * np.log1p(-1 / (4 * l[nz] ** 2))
        out[nz] = 2 * l[nz] * (X - Y_of_l(l[nz]))
    return out
def A0():
    """int_0^inf log(8 sqrt(u) Pi0_0(u)) du, u = v^4 near 0."""
    tot = 0.0
    # [0, 1]: u = v^4
    ed = np.array([0, 0.25, 0.5, 0.75, 1.0])
    v = np.concatenate([0.5 * (b - a) * xg + 0.5 * (a + b) for a, b in zip(ed[:-1], ed[1:])])
    W = np.concatenate([0.5 * (b - a) * wg for a, b in zip(ed[:-1], ed[1:])])
    u = v ** 4
    P0 = ff.Pi0(0, u)[0]
    tot += (np.log(8 * np.sqrt(u) * P0) * 4 * v ** 3) @ W
    ed = [1.0]
    while ed[-1] < 1e8:
        ed.append(ed[-1] * 2)
    ed = np.array(ed)
    u = np.concatenate([0.5 * (b - a) * xg + 0.5 * (a + b) for a, b in zip(ed[:-1], ed[1:])])
    W = np.concatenate([0.5 * (b - a) * wg for a, b in zip(ed[:-1], ed[1:])])
    th = ff.pi_th(0, u)[0]
    eps = -(2 / np.pi) * np.arctan(2 * ff.M0 / np.sqrt(u)) + 8 * np.sqrt(u) * th
    f = np.log1p(eps)
    b = f[-1] * u[-1] ** 1.5
    tot += f @ W + 2 * b * ed[-1] ** -0.5
    return tot




def delta0(u):
    """log(8 sqrt(u) Pi0_0(u)) without cancellation (u > 0)."""
    th = ff.pi_th(0, u)[0]
    sq = np.sqrt(u)
    eps = -(2 / np.pi) * np.arctan(2 * ff.M0 / sq) + 8 * sq * th
    return np.log1p(eps)


def int_delta0(ua):
    """int_0^ua delta0(u) du  (u = v^4 on [0,1], geometric panels above)."""
    ed = np.array([0, 0.25, 0.5, 0.75, 1.0])
    v = np.concatenate([0.5 * (b - a) * xg + 0.5 * (a + b) for a, b in zip(ed[:-1], ed[1:])])
    W = np.concatenate([0.5 * (b - a) * wg for a, b in zip(ed[:-1], ed[1:])])
    tot = (delta0(v ** 4) * 4 * v ** 3) @ W
    ed = [1.0]
    while ed[-1] * 1.5 < ua:
        ed.append(ed[-1] * 1.5)
    ed.append(ua)
    ed = np.array(ed)
    u = np.concatenate([0.5 * (b - a) * xg + 0.5 * (a + b) for a, b in zip(ed[:-1], ed[1:])])
    W = np.concatenate([0.5 * (b - a) * wg for a, b in zip(ed[:-1], ed[1:])])
    return tot + delta0(u) @ W


def R_of_J(J, r):
    """(2J+1)[log(P/Pi0) - C_J],  C_J = log(8 sqrt(u) Pi_vac(0,J)) = (1/2) log(1 - 1/(4 l^2)) + Y(l)."""
    J = np.asarray(J, float)
    K = 1.0 / r ** 2
    l = J + 0.5
    u = J * (J + 1) * K
    P0, P1, P2 = ff.Pi0(0, u)
    corr = -(K / 3) * u * P2 / P0
    C = 0.5 * np.log1p(-1 / (4 * l * l)) + Y_of_l(l)
    return (2 * J + 1) * (np.log1p(corr) - C)


def T_of_r(r, ts=100.0, J2fac=64):
    """T(r) = S0(r) - A0 r^2 - log(r)/6, never forming A0 itself:
    sum_{J>Js}(2J+1) delta0(u_J) = r^2 int_{ua}^inf delta0 du + Euler-Maclaurin boundary terms (midpoint form)."""
    Js = int(ts * r)
    J = np.arange(0, Js + 1, dtype=float)
    head = F(J, r).sum()
    a = Js + 0.5
    K = 1.0 / r ** 2
    ua = a * (a + 1) * K
    head -= int_delta0(ua) / K

    def G(x):
        return (2 * x + 1) * delta0(x * (x + 1) * K)
    h = 1.0
    xs = np.array([a - 2 * h, a - h, a + h, a + 2 * h])
    g = G(xs)
    Gp = (8 * (g[2] - g[1]) - (g[3] - g[0])) / (12 * h)
    Gppp = (g[3] - 2 * g[2] + 2 * g[1] - g[0]) / (2 * h ** 3)
    bdy = Gp / 24 - 7 * Gppp / 5760   # midpoint E-M: sum_{k>=0} f(a+1/2+k) = int_a f + f'(a)/24 - 7 f'''(a)/5760
    J2 = np.arange(Js + 1, J2fac * Js + 1, dtype=float)
    Rv = R_of_J(J2, r)
    rest = Rv.sum() + Rv[-1] * J2[-1] / 2.0   # sum_{J>N} c/J^3 ~ c/(2N^2) = R_N N/2
    return head + bdy + rest - np.log(r) / 6


if __name__ == '__main__':
    for r in [float(x) for x in sys.argv[1:]]:
        for ts in [60.0, 100.0, 160.0]:
            print(r, ts, repr(T_of_r(r, ts)), flush=True)
