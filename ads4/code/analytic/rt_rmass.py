"""
rt_rmass.py: certified unit-cell masses of the realised residual r = f - psi2 on [lo, hi] (cells of width 1), for the turning-point
band of Route T' (ANALYTIC_SLIVER §7).  f from CONTOUR's certified line panels (cb_common: Taylor polynomial of f about the panel centre,
truncation by Cauchy on a disc + far-RBF bound), psi2 subtracted exactly (Taylor shift), int |r| <= sum over subcells of
(poly_abs_bounds int of the r-polynomial) + 2 h e_f.  Beyond 1029: CONTOUR's tail(W, 0, T) (r = f for nu > Lam).
Output: cache/rt_rmass.json  [(T, T+1, mass), ...]
"""
import os, sys, json
sys.path.insert(0, '/home/user/MAJOR/ads4/code/contour')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cb_common as K, cb_lines as CL
from flint import arb, acb, acb_poly, ctx, fmpq_poly, fmpq
from an_common import nf02, LAM2, CACHE
ctx.prec = 128
W = K.load()
LAM = arb(LAM2).sqrt()
p2, _ = nf02()
Q2 = fmpq_poly([fmpq(c.numerator, c.denominator) for c in p2])


def psi_poly_at(t0, N):
    """Taylor coefficients in s of psi2((t0 + s)/Lam) (t0 integer + 1/2): exact shift in p then scaling by 1/Lam (ball)."""
    # psi2(p), p = (t0 + s)/Lam: coefficients of s^k = [shifted poly in p]_k / Lam^k
    p0 = acb(arb(t0)/LAM)
    # evaluate via Horner on acb_poly
    P = acb_poly([acb(arb(c.numerator)/c.denominator) for c in p2])
    lin = acb_poly([p0, acb(1/LAM)])
    return P(lin)


def cell_mass(t0, h=0.25, rho=1.0, N=40, Q=5):
    z0 = acb(t0)
    M1, F1 = K.G_disc_sup(W, z0, rho, 64)
    M2, F2 = K.G_disc_sup(W, z0.conjugate(), rho, 64)
    x = arb(h)/rho
    ef = (M1 + M2)*x**N/(1 - x) + F1 + F2
    pf = acb_poly(K.f_line_coeffs(W, z0, acb(1), N))
    pr = pf - psi_poly_at(t0, N)
    hs = arb(h)/Q
    tot = arb(0)
    for j in range(Q):
        cj = -arb(h) + (2*j + 1)*hs
        sr = pr(acb_poly([acb(cj), 1]))
        sup, intr = K.poly_abs_bounds(sr, hs)
        tot += intr + 2*hs*ef
    assert ef < arb('1e-14'), ef
    return float(tot.upper())


if __name__ == '__main__':
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    fn = os.path.join(CACHE, 'rt_rmass.json')
    out = json.load(open(fn)) if os.path.exists(fn) else []
    done = {r[0] for r in out}
    for T in range(lo, hi):
        if T in done:
            continue
        if T >= 1029:
            intf, _ = CL.tail(W, 0.0, T)
            m = float(intf.upper())
        else:
            m = cell_mass(T + 0.25) + cell_mass(T + 0.75)
        out.append([T, T + 1, m])
        print(T, m, flush=True)
        json.dump(sorted(out), open(fn, 'w'))
