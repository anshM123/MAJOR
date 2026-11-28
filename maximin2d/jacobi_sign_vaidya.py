"""Test of Proposition J (local crossing positivity from NEC via the cooperative Jacobi system)
in smooth planar Vaidya-AdS3:  ds^2 = z^-2(-(1-m(v)z^2)dv^2 - 2dvdz + dx^2),  NEC <=> m'(v) >= 0.

D(P,Q) = renormalized length of the (unique, here) spacelike geodesic between boundary points
P=(v1,x1), Q=(v2,x2) (cut-off z=EPS, D = s + 2 log EPS).  Prop. J predicts, for x1 < x2 and
null boundary directions  '+' = future-right (dv,dx)=(1,1), '-' = future-left (1,-1):
   D_{++} >= 0,  D_{--} >= 0  (same-type null directions at the two ends)
   D_{+-} <= 0,  D_{-+} <= 0  (opposite types; = 0 exactly in vacuum: D = log du + log dv)
and hence nu-density d_s d_t D >= 0 along any spacelike boundary slice.
Each mixed derivative is a 4-point finite difference.
usage: python3 -I jacobi_sign_vaidya.py
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

EPS = 2e-4


def make_mass(m0, dm, v0, w):
    def m(v):
        return m0 + dm * 0.5 * (1 + np.tanh((v - v0) / w))

    def mp(v):
        return dm * 0.5 / (w * np.cosh((v - v0) / w) ** 2)
    return m, mp


def rhs(s, y, m, mp):
    v, z, x, pv, pz, px = y
    f = 1 - m(v) * z * z
    dv = -z * z * pz
    dz = z * z * (f * pz - pv)
    dx = z * z * px
    dpv = 0.5 * z ** 4 * mp(v) * pz * pz
    # dH/dz, H = 1/2 z^2(-2 pv pz + f pz^2 + px^2), df/dz = -2 m z
    dpz = -(z * (-2 * pv * pz + f * pz * pz + px * px) + 0.5 * z * z * (-2 * m(v) * z) * pz * pz)
    return [dv, dz, dx, dpv, dpz, 0.0]


def shoot(v1, x1, E, p, m, mp, smax=60.0):
    z = EPS
    f = 1 - m(v1) * z * z
    disc = E * E - f * (p * p - 1 / z ** 2)
    pz = (E + np.sqrt(disc)) / f
    y0 = [v1, z, x1, E, pz, p]

    def hit(s, y, *a):
        return y[1] - EPS
    hit.terminal = True
    hit.direction = -1
    sol = solve_ivp(rhs, (0, smax), y0, args=(m, mp), method='DOP853', rtol=1e-12, atol=1e-14,
                    events=hit, first_step=1e-6)
    if sol.t_events[0].size == 0:
        return None
    ye = sol.y_events[0][0]
    return ye[0], ye[2], sol.t_events[0][0]


def D(v1, x1, v2, x2, m, mp, guess=None):
    dx = x2 - x1
    dt = v2 - v1
    # exact pure-AdS momenta: p_x = 2dx/(dx^2-dt^2), p_v = p_t = -2dt/(dx^2-dt^2)
    g0 = guess if guess is not None else np.array([-2 * dt, 2 * dx]) / (dx * dx - dt * dt)

    def F(q):
        r = shoot(v1, x1, q[0], q[1], m, mp)
        if r is None:
            return [1e3, 1e3]
        return [r[0] - v2, r[1] - x2]
    q, info, ier, msg = fsolve(F, g0, full_output=True, xtol=1e-13)
    r = shoot(v1, x1, q[0], q[1], m, mp)
    if r is None:
        return np.nan, q, np.inf
    res = max(abs(r[0] - v2), abs(r[1] - x2))
    return r[2] + 2 * np.log(EPS), q, res


def guess_by_continuation(P, Q, params, nsteps=8):
    m0, dm, v0, w = params
    q = None
    for lam in np.linspace(0, 1, nsteps + 1):
        mm, mmp = make_mass(lam * m0, lam * dm, v0, w)
        d, q, res = D(P[0], P[1], Q[0], Q[1], mm, mmp, q)
        assert res < 1e-9, (lam, res)
    return q


def mixed(P, Q, dP, dQ, h, params):
    m, mp = make_mass(*params)
    vals = {}
    q = guess_by_continuation(P, Q, params)
    for a in (1, -1):
        for b in (1, -1):
            d, q, res = D(P[0] + a * h * dP[0], P[1] + a * h * dP[1],
                          Q[0] + b * h * dQ[0], Q[1] + b * h * dQ[1], m, mp, q)
            assert res < 1e-9, res
            vals[(a, b)] = d
    return (vals[(1, 1)] - vals[(1, -1)] - vals[(-1, 1)] + vals[(-1, -1)]) / (4 * h * h)


if __name__ == '__main__':
    plus, minus = (1.0, 1.0), (1.0, -1.0)
    spat = lambda tp: (tp, 1.0)   # spacelike boundary tangent (dv,dx) = (t',1), |t'|<1
    h = 0.02
    cases = [("vacuum planar BTZ m=1", (1.0, 0.0, 0.0, 0.3)),
             ("NEC shell dm=+1.5", (1.0, 1.5, 0.0, 0.3)),
             ("NEC shell dm=+4", (0.5, 4.0, 0.2, 0.25)),
             ("NEC-violating dm=-0.8", (1.0, -0.8, 0.0, 0.3)),
             ("NEC-violating dm=-0.95", (1.0, -0.95, 0.2, 0.25))]
    configs = [((0.0, 0.0), (0.4, 1.5)), ((-0.3, 0.0), (0.6, 1.2)), ((0.5, 0.0), (-0.2, 1.3)),
               ((0.2, 0.0), (0.2, 2.0))]
    for name, (m0, dm, v0, w) in cases:
        prm = (m0, dm, v0, w)
        print("==", name)
        for P, Q in configs:
            Dpp = mixed(P, Q, plus, plus, h, prm)
            Dmm = mixed(P, Q, minus, minus, h, prm)
            Dpm = mixed(P, Q, plus, minus, h, prm)
            Dmp = mixed(P, Q, minus, plus, h, prm)
            Dss = mixed(P, Q, spat(0.5), spat(-0.3), h, prm)
            print(f"  P={P} Q={Q}:  D++={Dpp:+.3e} D--={Dmm:+.3e} D+-={Dpm:+.3e} D-+={Dmp:+.3e}"
                  f"  nu-density(spacelike)={Dss:+.3e}")
