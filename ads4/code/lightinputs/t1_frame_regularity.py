"""t1: polarisation factors in the CHMRS Regge frame add no powers of r  (LIGHT_INPUTS.md Sec. 2, Lemma F).

Frame (CHMRS (C.2), SPINNING_E2E.md 1.1): (x1,x2,x3,x4) = (-y, x/x^2, y, -x/x^2), x,y future timelike, r=|x||y|,
eta = -xhat.yhat.  We fix |x| = 1, so x2 = -xhat, x4 = +xhat, x1 = -r yhat, x3 = r yhat.
Polarisations: at 1,3 the frame circle z(th) = yhat + cos th e_par + sin th e_perp (r-independent);
at 2,4 the inversion-referred circle I(xhat) z'(th) (r-independent).
We evaluate every weight-0 conformal building block (normalised V_{i,jk}, H_ij) at r -> 0 and fit the power of r.
Claim [DV]: every block tends to a finite limit (exponent 0, or a positive integer if it vanishes), i.e. the
stripped tensor s13^D s24^D <...> = u^{-D} sum_I tau_I(z;r,eta) f_I(u,v) has tau_I analytic at r = 0: the Regge
growth of every frame mode is that of the scalar-like functions u^{-D} f_I, i.e. r^{1-J0} with J0 the
Regge intercept as defined by KKSDZ (1904.05905, eq. J0definition: fixed generic polarisations, finite in the limit).
"""
import numpy as np
from itertools import permutations

eta_metric = np.diag([-1.0, 1.0, 1.0])
dot = lambda a, b: a @ eta_metric @ b

def frame(eta):
    a = np.arccosh(eta)
    yh = np.array([1.0, 0.0, 0.0])
    xh = np.array([np.cosh(a), np.sinh(a), 0.0])
    # unit tangent at yhat pointing to xhat, and the third unit vector
    epar = (xh - eta * yh) / np.sqrt(eta**2 - 1)
    eperp = np.array([0.0, 0.0, 1.0])
    # tangent at xhat pointing to yhat
    fpar = (yh - eta * xh) / np.sqrt(eta**2 - 1)
    return xh, yh, epar, eperp, fpar

def config(r, eta, th):
    xh, yh, epar, eperp, fpar = frame(eta)
    X = [-r * yh, -xh, r * yh, xh]
    z13 = lambda t: yh + np.cos(t) * epar + np.sin(t) * eperp
    I = np.eye(3) + 2 * np.outer(xh, xh @ eta_metric)      # inversion tensor I(xhat)^mu_nu, xhat^2 = -1
    z24 = lambda t: I @ (xh + np.cos(t) * fpar + np.sin(t) * eperp)
    Z = [z13(th[0]), z24(th[1]), z13(th[2]), z24(th[3])]
    for z in Z:
        assert abs(dot(z, z)) < 1e-12
    return X, Z

def blocks(X, Z):
    PP = lambda i, j: -dot(X[i] - X[j], X[i] - X[j]) / 2
    ZP = lambda i, j: dot(Z[i], X[j] - X[i])
    ZZ = lambda i, j: dot(Z[i], Z[j])
    out = {}
    for i in range(4):
        for j, k in permutations([m for m in range(4) if m != i], 2):
            if j < k:
                V = (ZP(i, j) * PP(i, k) - ZP(i, k) * PP(i, j)) / PP(j, k)
                out[f"V{i+1},{j+1}{k+1}"] = V / np.sqrt(abs(PP(i, j) * PP(i, k) / PP(j, k)))
    for i in range(4):
        for j in range(i + 1, 4):
            H = -2 * (ZZ(i, j) * PP(i, j) - ZP(i, j) * ZP(j, i))
            out[f"H{i+1}{j+1}"] = H / PP(i, j)
    return out

def main():
    rng = np.random.default_rng(1)
    rs = np.array([1e-2, 1e-3, 1e-4, 1e-5])
    worst = 0.0
    print("t1: weight-0 building blocks in the CHMRS frame as r -> 0 (fixed eta, fixed frame polarisations)")
    for eta in [1.3, 2.0, 5.0]:
        th = rng.uniform(0, 2 * np.pi, 4)
        vals = np.array([[v for v in blocks(*config(r, eta, th)).values()] for r in rs])
        names = list(blocks(*config(rs[0], eta, th)).keys())
        lim = vals[-1]
        # power fit on |B(r) - B(0)| tells analyticity order; on |B(r)| tells the leading power
        print(f"eta = {eta}, theta = {np.round(th,3)}")
        for n, col in zip(names, vals.T):
            p = np.polyfit(np.log(rs), np.log(np.abs(col) + 1e-300), 1)[0]
            dif = np.abs(col[:-1] - col[-1])
            q = np.polyfit(np.log(rs[:-1]), np.log(dif + 1e-300), 1)[0] if dif.max() > 1e-13 else np.inf
            worst = min(worst, p)
            print(f"   {n:9s} B(r=1e-5) = {col[-1]: .6f}   leading power {p: .3f}   B(r)-B(0) ~ r^{q:.2f}")
    print(f"most negative leading power over all blocks/eta: {worst:.3e}  (0 = finite limit; no negative power)")

if __name__ == "__main__":
    main()
