"""SO(2)-symmetric extremal surfaces in Vaidya-AdS4 (Poincare, ingoing EF).

  ds^2 = z^-2 ( -f dv^2 - 2 dv dz + drho^2 + rho^2 dphi^2 ),  f = 1 - m(v) z^3,
  m(v) = M/2 (1 + tanh(v/dlt))     (thick shell; m' >= 0  <=>  NEC)
  dlt -> 0 is the thin shell.

A rotationally symmetric surface is a curve x(s)=(v,z,rho) in the orbit space;
Area/(2 pi) = G-length with G = W h, W = rho^2/z^4, h = -f dv^2 - 2 dv dz + drho^2.
We integrate G-geodesics in the h-unit-speed parameter sigma with the rescaled
momenta pi = h(u,.), u = dx/dsigma:

  u^v = -pi_z, u^z = -pi_v + f pi_z, u^rho = pi_rho,
  dpi_a/dsigma = -(1/2) d_a f pi_z^2 + (q/2) d_a lnW - (1/2) pi_a (u . d lnW),
  dA/dsigma = rho/z^2,     q = h^{-1}(pi,pi) = 1 (conserved).

Renormalised area (per 2 pi): A_ren = lim (A_eps - sum_anchors rho_p/eps).
Near an anchor rho = rho_p - z^2/(2 rho_p) + al z^3, and the integrand is
rho_p/z^2 + z (m rho_p/2 - 2 al) + O(z^2); this O(eps^2) remainder is added at both ends.
"""
import numpy as np
from scipy.integrate import solve_ivp


class Vaidya:
    def __init__(self, M=1.0, dlt=0.1, NEC=True):
        self.M = M; self.dlt = dlt; self.sign = 1.0 if NEC else -1.0

    def m(self, v):
        # NEC=False: negative-energy shell, m(v) = -M/2 (1+tanh) (m' < 0 violates the NEC)
        return self.sign * 0.5 * self.M * (1.0 + np.tanh(v / self.dlt))

    def mp(self, v):
        c = 1.0 / np.cosh(np.clip(v / self.dlt, -350, 350))
        return self.sign * 0.5 * self.M / self.dlt * c * c

    # ---------------- vector field (works for scalars and numpy arrays) -------------
    def rhs(self, s, y):
        v, z, r, pv, pz, pr, A = y
        mv = self.m(v); f = 1.0 - mv * z**3
        dfv = -self.mp(v) * z**3; dfz = -3.0 * mv * z**2
        uv = -pz; uz = -pv + f * pz; ur = pr
        q = -2.0 * pv * pz + f * pz * pz + pr * pr
        gz = -4.0 / z; gr = 2.0 / r
        ud = uz * gz + ur * gr
        dpv = -0.5 * dfv * pz * pz - 0.5 * pv * ud
        dpz = -0.5 * dfz * pz * pz + 0.5 * q * gz - 0.5 * pz * ud
        dpr = 0.5 * q * gr - 0.5 * pr * ud
        return [uv, uz, ur, dpv, dpz, dpr, r / z**2]

    def f(self, v, z):
        return 1.0 - self.m(v) * z**3

    def pi_from_u(self, v, z, u):
        uv, uz, ur = u
        f = self.f(v, z)
        n = -f * uv * uv - 2 * uv * uz + ur * ur
        s = 1.0 / np.sqrt(n)
        uv, uz, ur = uv * s, uz * s, ur * s
        return -f * uv - uz, -uv, ur

    # ---------------- initial data ----------------
    def cap_init(self, vs, zs, d0=1e-4):
        f = self.f(vs, zs)
        v = vs + d0 * d0 / (2 * zs); z = zs - f * d0 * d0 / (2 * zs); r = d0
        pv, pz, pr = self.pi_from_u(v, z, (d0 / zs, -f * d0 / zs, 1.0))
        return np.array([v, z, r, pv, pz, pr, d0 * d0 / (2 * zs * zs)])

    def edge_init(self, r1, t1, al, be, zs=1e-3):
        v = t1 - zs + be * zs**3; z = zs; r = r1 - zs * zs / (2 * r1) + al * zs**3
        pv, pz, pr = self.pi_from_u(v, z, (-1 + 3 * be * zs * zs, 1.0, -zs / r1 + 3 * al * zs * zs))
        # renormalised contribution of [0,zs]: integrand = r1/z^2 + z (m r1/2 - 2 al) + O(z^2)
        A0 = -r1 / zs + zs * zs * (self.m(t1) * r1 / 4 - al)
        return np.array([v, z, r, pv, pz, pr, A0])

    # ---------------- precise integration ----------------
    def integrate(self, y0, zend=1e-3, smax=60.0, zmax=25.0, rmin=2e-4, rtol=1e-11, atol=1e-12):
        def ev_bdy(s, y):
            return y[1] - zend
        ev_bdy.terminal = True; ev_bdy.direction = -1
        def ev_deep(s, y):
            return y[1] - zmax
        ev_deep.terminal = True
        def ev_axis(s, y):
            return y[2] - rmin
        ev_axis.terminal = True; ev_axis.direction = -1
        sol = solve_ivp(self.rhs, (0, smax), y0, method='DOP853', rtol=rtol, atol=atol,
                        events=(ev_bdy, ev_deep, ev_axis), dense_output=False)
        y = sol.y[:, -1]
        if sol.status == 1 and len(sol.t_events[0]):
            r = y[2]; z = y[1]
            uv, uz, ur = self.rhs(0, y)[:3]
            drdz, dvdz = ur / uz, uv / uz
            r0 = r + z * z / (2 * r)
            al2 = (drdz + z / r0) / (3 * z * z)
            be2 = (dvdz + 1) / (3 * z * z)
            r_anchor = r + z * z / (2 * r0) - al2 * z**3
            t_anchor = y[0] + z - be2 * z**3
            A = y[6] - r_anchor / z + z * z * (self.m(t_anchor) * r_anchor / 4 - al2)
            return dict(status='bdy', r=r_anchor, t=t_anchor, A=A, y=y, s=sol.t[-1])
        if sol.status == 1 and len(sol.t_events[1]):
            return dict(status='deep', y=y)
        if sol.status == 1 and len(sol.t_events[2]):
            return dict(status='axis', y=y)
        return dict(status='other', y=y)

    def cap(self, vs, zs, **kw):
        return self.integrate(self.cap_init(vs, zs), **kw)

    def tube(self, r1, t1, al, be, zs=1e-3, **kw):
        return self.integrate(self.edge_init(r1, t1, al, be, zs=zs), zend=zs, **kw)

    # ---------------- coarse batch integration (fixed-step RK4) ----------------
    def batch(self, Y0, h=0.004, nmax=40000, zstop=0.03, zmax=25.0, rmin=0.005):
        """Y0: (7,N).  Returns status array and landing (r,t) estimates (coarse)."""
        Y = Y0.copy(); N = Y.shape[1]
        alive = np.ones(N, bool); status = np.zeros(N, int)  # 0 run,1 bdy,2 deep,3 axis,4 timeout
        out = np.full((3, N), np.nan)
        started = np.zeros(N, bool)  # must first leave the near-boundary layer
        for it in range(nmax):
            idx = np.nonzero(alive)[0]
            if len(idx) == 0:
                break
            y = Y[:, idx]
            k1 = np.array(self.rhs(0, y)); k2 = np.array(self.rhs(0, y + 0.5 * h * k1))
            k3 = np.array(self.rhs(0, y + 0.5 * h * k2)); k4 = np.array(self.rhs(0, y + h * k3))
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            Y[:, idx] = y
            z = y[1]; r = y[2]
            started[idx] |= z > 2 * zstop
            bad = ~np.all(np.isfinite(y), axis=0)
            land = started[idx] & (z < zstop) & (y[4] * 0 == 0)
            deep = z > zmax
            ax = r < rmin
            for code, msk in ((1, land), (2, deep | bad), (3, ax & ~land)):
                sel = idx[msk & alive[idx]]
                status[sel] = code; alive[sel] = False
            sel = idx[land]
            out[0, sel] = y[2, land]; out[1, sel] = y[0, land] + y[1, land]; out[2, sel] = y[6, land]
        status[alive] = 4
        return status, out
