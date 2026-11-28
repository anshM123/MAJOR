"""All extremal pieces (caps and tubes) for concentric annuli in Vaidya-AdS4.

caps : shoot from the axis with data (v*, z*); coarse batch over a grid, find
       every grid triangle whose image contains the target (R,t), refine by
       Newton with the precise integrator.  c(R,t) = min area over branches.
tubes: shoot from edge (rho1,t1) with the free z^3 data (alpha,beta); same
       triangle search for the target (rho2,t2).  Every branch is returned.
Areas are renormalised per 2 pi (A_ren = A - sum rho_p/eps).
"""
import numpy as np
from vaidya4 import Vaidya


def batch_adaptive(V, Y0, hfac=0.02, nmax=8000, zstop=0.0008, zmax=6.0, rmin=3e-4, smax=30.0):
    Y = Y0.copy(); N = Y.shape[1]
    alive = np.ones(N, bool); status = np.zeros(N, int)
    out = np.full((3, N), np.nan); s = np.zeros(N)
    zmax_seen = Y[1].copy(); rmax_seen = Y[2].copy()
    for it in range(nmax):
        idx = np.nonzero(alive)[0]
        if len(idx) == 0:
            break
        y = Y[:, idx]
        h = hfac * np.minimum(np.minimum(y[1], np.abs(y[2])), 1.0)
        k1 = np.array(V.rhs(0, y)); k2 = np.array(V.rhs(0, y + 0.5 * h * k1))
        k3 = np.array(V.rhs(0, y + 0.5 * h * k2)); k4 = np.array(V.rhs(0, y + h * k3))
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        Y[:, idx] = y; s[idx] += h
        zmax_seen[idx] = np.maximum(zmax_seen[idx], y[1]); rmax_seen[idx] = np.maximum(rmax_seen[idx], y[2])
        z, r = y[1], y[2]
        bad = ~np.all(np.isfinite(y), axis=0)
        land = (z < zstop) & (zmax_seen[idx] > 2.5 * zstop)
        deep = (z > zmax) | bad | (s[idx] > smax)
        ax = (r < rmin) & (rmax_seen[idx] > 4 * rmin) & ~land
        for code, msk in ((1, land), (2, deep & ~land), (3, ax & ~deep)):
            sel = idx[msk]
            sel = sel[alive[sel]]
            status[sel] = code; alive[sel] = False
        sel = idx[land]
        rr = y[2, land]; zz = y[1, land]
        out[0, sel] = rr + zz * zz / (2 * rr); out[1, sel] = y[0, land] + zz; out[2, sel] = y[6, land]
    status[alive] = 4
    return status, out


def tri_candidates(P1, P2, ok, X, Y, target, maxdiam=0.3, bigcells=None, bigdiam=6.0):
    """Grid (n1 x n2) params P1,P2, images X,Y (landing r,t), mask ok.  Return
    param points (barycentric interpolation) of triangles whose image contains target."""
    n1, n2 = X.shape
    cands = []
    tx, ty = target
    for (a, b, c) in (((0, 0), (1, 0), (0, 1)), ((1, 1), (0, 1), (1, 0))):
        sl = lambda o: (slice(o[0], n1 - 1 + o[0]), slice(o[1], n2 - 1 + o[1]))
        xa, ya, xb, yb, xc, yc = X[sl(a)], Y[sl(a)], X[sl(b)], Y[sl(b)], X[sl(c)], Y[sl(c)]
        m = ok[sl(a)] & ok[sl(b)] & ok[sl(c)]
        det = (xb - xa) * (yc - ya) - (xc - xa) * (yb - ya)
        with np.errstate(all='ignore'):
            l1 = ((tx - xa) * (yc - ya) - (xc - xa) * (ty - ya)) / det
            l2 = ((xb - xa) * (ty - ya) - (tx - xa) * (yb - ya)) / det
        diam = np.maximum.reduce([np.hypot(xb - xa, yb - ya), np.hypot(xc - xa, yc - ya), np.hypot(xc - xb, yc - yb)])
        inside0 = m & (l1 >= -1e-9) & (l2 >= -1e-9) & (l1 + l2 <= 1 + 1e-9) & np.isfinite(det) & (det != 0)
        if bigcells is not None:
            for i, j in zip(*np.nonzero(inside0 & (diam >= maxdiam) & (diam < bigdiam))):
                bigcells.add((i, j))
        inside = inside0 & (diam < maxdiam)
        for i, j in zip(*np.nonzero(inside)):
            pa = np.array([P1[i + a[0], j + a[1]], P2[i + a[0], j + a[1]]])
            pb = np.array([P1[i + b[0], j + b[1]], P2[i + b[0], j + b[1]]])
            pc = np.array([P1[i + c[0], j + c[1]], P2[i + c[0], j + c[1]]])
            L1, L2 = l1[i, j], l2[i, j]
            cands.append((1 - L1 - L2) * pa + L1 * pb + L2 * pc)
    return cands


def newton2(fun, x0, target, tol=1e-9, maxit=30, h=1e-6):
    x = np.array(x0, float)
    for it in range(maxit):
        r0 = fun(x)
        if r0 is None:
            return None
        F0 = np.array(r0[:2]) - target
        if np.max(np.abs(F0)) < tol:
            return x, r0
        J = np.zeros((2, 2))
        for k in range(2):
            dx = np.zeros(2); dx[k] = h * max(1.0, abs(x[k]))
            r1 = fun(x + dx)
            if r1 is None:
                return None
            J[:, k] = (np.array(r1[:2]) - np.array(r0[:2])) / dx[k]
        try:
            step = np.linalg.solve(J, -F0)
        except np.linalg.LinAlgError:
            return None
        lam = 1.0
        for _ in range(12):
            r2 = fun(x + lam * step)
            if r2 is not None and np.max(np.abs(np.array(r2[:2]) - target)) < np.max(np.abs(F0)):
                break
            lam *= 0.5
        else:
            # no decrease: accept if already at the integration noise floor
            if np.max(np.abs(F0)) < 1e-8:
                return x, r0
            return None
        x = x + lam * step
    r0 = fun(x)
    if r0 is not None and np.max(np.abs(np.array(r0[:2]) - target)) < 1e-8:
        return x, r0
    return None


def cluster(cands, rel=1e-2):
    out = []
    for c in cands:
        if not any(np.max(np.abs(c - d) / (1e-3 + np.abs(d))) < rel for d in out):
            out.append(c)
    return out


def continuation(fun, x0, y0, y1, h0=0.05, hmin=1e-4):
    """Follow a solution of fun(x)[:2] = y along the straight path y0 -> y1."""
    x = np.array(x0, float); s = 0.0; h = h0; res = None
    y0 = np.array(y0, float); y1 = np.array(y1, float)
    while s < 1.0 - 1e-14:
        sn = min(1.0, s + h)
        r = newton2(fun, x, y0 + sn * (y1 - y0), maxit=12)
        if r is None:
            h *= 0.5
            if h < hmin:
                return None
            continue
        x, res = r; s = sn; h = min(2 * h, 0.2)
    return x, res


def dedupe(sols, tol=1e-6):
    out = []
    for x, r in sols:
        if not any(abs(r[2] - r2[2]) < tol and np.max(np.abs(x - x2)) < 1e-4 * max(1, np.max(np.abs(x))) for x2, r2 in out):
            out.append((x, r))
    return out


def find_preimages(V, initfun, P1, P2, ok, X, Y, target, depth=2, nsub=9):
    """Triangle search with local batch refinement of cells whose image is large."""
    big = set()
    cands = tri_candidates(P1, P2, ok, X, Y, target, bigcells=big)
    if depth == 0 or not big:
        return cands
    big = sorted(big)[:400]
    subs = []
    for (i, j) in big:
        a0, a1 = P1[i, j], P1[i + 1, j]; b0, b1 = P2[i, j], P2[i, j + 1]
        A, B = np.meshgrid(np.linspace(a0, a1, nsub), np.linspace(b0, b1, nsub), indexing='ij')
        subs.append((A, B))
    Y0 = np.array([initfun(a, b) for A, B in subs for a, b in zip(A.ravel(), B.ravel())]).T
    st, out = batch_adaptive(V, Y0)
    n = nsub * nsub
    for k, (A, B) in enumerate(subs):
        sl = slice(k * n, (k + 1) * n)
        okk = (st[sl] == 1).reshape(A.shape)
        cands += find_preimages(V, initfun, A, B, okk, out[0, sl].reshape(A.shape), out[1, sl].reshape(A.shape),
                                target, depth - 1, nsub)
    return cands


class Pieces:
    def __init__(self, V, cap_grid=(200, 200), tube_grid=(321, 321), zs=2e-4):
        self.V = V; self.zs = zs; self.cap_grid = cap_grid; self.tube_grid = tube_grid
        self._capscan = None; self._tubescan = {}

    # ---------- caps ----------
    def cap_fun(self, x):
        vs, zs = x
        if zs <= 1e-3:
            return None
        r = self.V.cap(vs, zs)
        if r['status'] != 'bdy':
            return None
        return (r['r'], r['t'], r['A'])

    def cap_scan(self, vrange=(-7.0, 4.0), zrange=(0.02, 9.0)):
        n1, n2 = self.cap_grid
        vs = np.linspace(*vrange, n1)
        # dense near the late-time horizon z_h = M^(-1/3) where large discs have their tip
        zh = self.V.M ** (-1.0 / 3) if self.V.M > 0 else 1.0
        zz = np.unique(np.concatenate([np.geomspace(*zrange, n2),
                                       zh * (1 - np.geomspace(1e-10, 0.2, n2 // 2)),
                                       zh * (1 + np.geomspace(1e-6, 0.2, n2 // 4))]))
        P1, P2 = np.meshgrid(vs, zz, indexing='ij')
        Y0 = np.array([self.V.cap_init(a, b, d0=1e-3) for a, b in zip(P1.ravel(), P2.ravel())]).T
        st, out = batch_adaptive(self.V, Y0)
        ok = (st == 1).reshape(P1.shape)
        self._capscan = (P1, P2, ok, out[0].reshape(P1.shape), out[1].reshape(P1.shape))
        return self._capscan

    def caps(self, R, t, continue_from_vacuum=True):
        if self._capscan is None:
            self.cap_scan()
        P1, P2, ok, X, Y = self._capscan
        sols = []
        if continue_from_vacuum:
            ts = min(t, -8 * self.V.dlt - 0.05)
            r = continuation(self.cap_fun, (ts - R, R), (R, ts), (R, t))
            if r is not None:
                sols.append(r)
        init = lambda a, b: self.V.cap_init(a, b, d0=1e-3)
        for x0 in cluster(find_preimages(self.V, init, P1, P2, ok, X, Y, (R, t))):
            res = newton2(self.cap_fun, x0, np.array([R, t]))
            if res is not None:
                sols.append(res)
        return dedupe(sols)

    # ---------- tubes ----------
    def tube_fun(self, r1, t1):
        zs = self.zs
        def fun(x):
            # validity of the near-boundary expansion and no spurious immediate return
            if abs(x[0]) * zs**2 > 0.05 * r1**-1 or abs(x[1]) * zs**2 > 0.05 * r1**-1:
                return None
            res = self.V.tube(r1, t1, x[0], x[1], zs=zs)
            if res['status'] != 'bdy' or res['s'] < 20 * zs:
                return None
            return (res['r'], res['t'], res['A'])
        return fun

    def tube_scan(self, r1, t1, urange=13.0):
        key = (round(r1, 12), round(t1, 12))
        if key in self._tubescan:
            return self._tubescan[key]
        n1, n2 = self.tube_grid
        u = np.linspace(-urange, urange, n1); w = np.linspace(-urange, urange, n2)
        P1, P2 = np.meshgrid(np.sinh(u) / r1**3, np.sinh(w) / r1**3, indexing='ij')
        Y0 = np.array([self.V.edge_init(r1, t1, a, b, zs=0.002) for a, b in zip(P1.ravel(), P2.ravel())]).T
        st, out = batch_adaptive(self.V, Y0)
        ok = (st == 1).reshape(P1.shape)
        scan = (P1, P2, ok, out[0].reshape(P1.shape), out[1].reshape(P1.shape))
        self._tubescan[key] = scan
        return scan

    def tubes(self, r1, t1, r2, t2):
        P1, P2, ok, X, Y = self.tube_scan(r1, t1)
        sols = []
        init = lambda a, b: self.V.edge_init(r1, t1, a, b, zs=0.002)
        for x0 in cluster(find_preimages(self.V, init, P1, P2, ok, X, Y, (r2, t2))):
            res = newton2(self.tube_fun(r1, t1), x0, np.array([r2, t2]))
            if res is not None:
                sols.append(res)
        return dedupe(sols)
