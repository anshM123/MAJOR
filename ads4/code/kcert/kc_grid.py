"""kc_grid.py -- exact (ball-arithmetic) light constants of the GRID form of Omega^EE's k = 2 piece.

Grid functional (prl/data/README.md Sec. 4, 7, 8): Re sum_{nu_k <= Lam+10} dnu X_k Psi^_(EE,nu_k)/(h_- + h_+)(nu_k), X_k as printed
(decimal strings, read exactly).  Its light values:
  |g3|^2 (published "unshifted" reading):  sum_k dnu Re X_k p_k^6 / Lam   ->  dK0^grid = -2 [ that - int_0^1 psi2 p^6 dp ]
  Einstein (V3):  Re sum_k dnu chi_k [R_-(nu_k) D(nu_k + i) + R_+(nu_k) D(nu_k - i)],  chi_k = X_k/(R_- + R_+)(nu_k),  vs L_A/Lam.
Files: d_abel_X_{300,1000}_trunc.txt (dnu = 0.02 / 0.04), gap/d_abel_X_{2000,5000}_trunc.txt.gz (dnu = 0.05).
The sums are exact up to ball rounding (no quadrature error: the grid sum IS the functional).
Usage: python3 kc_grid.py [Dg ...]   -> cache/kc_grid.txt
"""
import sys, os, gzip, time
from fractions import Fraction as Fr
from flint import arb, acb, ctx, fmpq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kc_common import *  # noqa
from kc_dk0 import LA_over_Lam

ctx.prec = 128
FILES = {300: ('d_abel_X_300_trunc.txt', Fr(2, 100)), 1000: ('d_abel_X_1000_trunc.txt', Fr(4, 100)),
         2000: ('gap/d_abel_X_2000_trunc.txt.gz', Fr(5, 100)), 5000: ('gap/d_abel_X_5000_trunc.txt.gz', Fr(5, 100))}


def rows(fn):
    op = gzip.open if fn.endswith('.gz') else open
    with op(os.path.join(DATA, fn), 'rt') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            yield line.split()


def main(Dgs):
    out = []; res = {}

    def pr(s):
        print(s, flush=True); out.append(s)
    p2, p3 = nf02(); K0 = K0_exact()
    Rm_s, Rp_s = kernels()
    Rm, Rp = Rat(Rm_s), Rat(Rp_s)
    I = acb(0, 1)
    for Dg in Dgs:
        t0 = time.time()
        fn, dnu = FILES[Dg]
        L = Lam(Dg); L2 = arb(lam_fmpq_sq(Dg))
        s6 = arb(0); v3 = arb(0); s6s = arb(0); n = 0; numax = 0.0
        for r in rows(fn):
            nu = arb(r[0]); ReX = arb(r[1]); ImX = arb(r[2])
            k = Fr(r[0]) / dnu - Fr(1, 2)
            assert k.denominator == 1, r[0]          # nu_k = (k + 1/2) dnu exactly
            nu2 = nu * nu
            s6 += ReX * nu2 * nu2 * nu2
            z = acb(nu)
            Rmz, Rpz = Rm(z), Rp(z)
            chi = acb(ReX, ImX) / (Rmz + Rpz)
            s6s += (chi * (Rmz * (z + I) ** 6 + Rpz * (z - I) ** 6)).real
            v3 += (chi * (Rmz / ((z + I) ** 2 + arb(9) / 4) + Rpz / ((z - I) ** 2 + arb(9) / 4))).real
            n += 1; numax = float(r[0])
        assert numax <= float(L) + 10.0 + 1e-9
        d = arb(fq(dnu))
        J6 = d * s6 / (L2 * L2 * L2 * L)
        dK0 = -2 * (J6 - arb(fq(poly_moment(p2, 6))))
        V3 = d * v3
        dK0s = -2 * (d * s6s / (L2 * L2 * L2 * L) - arb(fq(poly_moment(p2, 6))))
        LA = LA_over_Lam(Dg, p2)
        pr('Dg = %d (%s, %d rows, nu <= %.3f): dK0^grid = %s ; Lam^2 dK0 = %s ; K0 + dK0 = %s' % (
            Dg, fn, n, numax, dK0.str(12, radius=True), (dK0 * L2).str(8, radius=True), (arb(fq(K0)) + dK0).str(12, radius=True)))
        pr('   V3^grid = %s vs L_A/Lam = %s (rel diff %s); 16 L_A = %s  [%.0fs]' % (
            V3.str(15, radius=True), (LA / L).str(15), ((V3 - LA / L) / (LA / L)).str(3, radius=True), (16 * LA).str(12), time.time() - t0))
        pr('   shifted (C-prime) reading: dK0s^grid = %s ; Lam^2 dK0s = %s' % (dK0s.str(8, radius=True), (dK0s * L2).str(6, radius=True)))
        res[Dg] = dict(dK0=dK0, dK0s=dK0s, V3=V3, LA=LA, Lam=L, K0=arb(fq(K0)))
    return out, res


if __name__ == '__main__':
    Dgs = [int(a) for a in sys.argv[1:]] or [300, 1000, 2000, 5000]
    o, _ = main(Dgs)
    open(os.path.join(HERE, 'cache', 'kc_grid.txt'), 'w').write('\n'.join(o) + '\n')
