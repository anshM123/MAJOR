r"""Certificate C4:   FF(b, beta) := sum_{n>=0} beta^{2n} G_n(b) > 0   on the region
        R = { 0 <= b <= TAU2 = 10.85,   0 <= beta <= b / MK },     MK = 6,
where F_M(beta) = FF(M beta, beta) for beta < pi (Lemma 7).  R contains every (b, beta) = (M beta, beta) with
M >= MK and M beta <= TAU2.   The slice beta = 0 is the flat functional fhat(b) = G_0(b).

For a b-cell [bc - h, bc + h] and x = beta^2 in [x0, x1] (x1 <= (TAU2/MK)^2 < pi^2) we use, for all b in the cell,
   FF(b, beta) >= P_c(x) - Mu_N(x1) * ( (283/120) h + (2591/420) tailK(bc + h) ) - (2591/420) * tailN(x1),
   P_c(x)   = sum_{n<=N} x^n G_n^{(K)}(bc),     G_n^{(K)}(bc) = sum_{k<=K} g_{n,k} bc^{2k}    (exact rationals g_{n,k}),
   Mu_N(x)  = sum_{n<=N} mu_n x^n   (mu_n = [y^{2n}] y/sin y >= 0),
   tailK(b) = (b^{2K+2}/(2K+2)!) / (1 - b^2/((2K+3)(2K+4))),
   tailN(x) = 2 (x/pi^2)^{N+1} / (1 - x/pi^2).
P_c is enclosed on [x0, x1] by interval Horner evaluation and by the centred (mean-value) form; x-intervals are bisected
adaptively, b-cells are bisected when the x-bisection does not terminate.  One output line per certified box.
"""
import sys
from fractions import Fraction as Fr
from math import factorial
from multiprocessing import Pool
from common import *
from expansion import g_coeffs, mu_coeffs

N = 14
K = 36
G = g_coeffs(N, K)
MU = mu_coeffs(N)
G_IV = [[frac_iv(G[n][k]) for k in range(K + 1)] for n in range(N + 1)]
MU_IV = [frac_iv(m) for m in MU]
NF = frac_iv(NORM_F)
NPF = frac_iv(NORM_PF)
PI2_LO = (iv.pi * iv.pi).a
XMAX_GLOBAL = (TAU2 / MK) ** 2
assert frac_iv(XMAX_GLOBAL).b < PI2_LO


def horner(coefs, x):
    v = iv.mpf(0)
    for a in reversed(coefs):
        v = v * x + a
    return v


def tailK(b_hi):
    b2 = frac_iv(b_hi) ** 2
    return ipow(b2, K + 1) / factorial(2 * K + 2) / (1 - b2 / ((2 * K + 3) * (2 * K + 4)))


def tailN(x1):
    r = x1 / iv.mpf(PI2_LO)
    return 2 * ipow(r, N + 1) / (1 - r)


def certify_cell(args):
    """returns (ok, lines, min_lower_bound) for the b-cell [b0, b1]."""
    b0, b1 = args
    bc, h = (b0 + b1) / 2, (b1 - b0) / 2
    y = frac_iv(bc) ** 2
    Gc = [horner(G_IV[n], y) for n in range(N + 1)]                 # enclosures of G_n^{(K)}(bc)
    dGc = [(n + 1) * Gc[n + 1] for n in range(N)]                    # coefficients of P_c'
    errB = NPF * frac_iv(h) + NF * tailK(b1)
    X = (b1 / MK) ** 2
    lines = []
    stack = [(Fr(0), X, 0)]
    vmin = None
    nx = 0
    while stack:
        x0, x1, depth = stack.pop()
        xi = iv.mpf([frac_iv(x0).a, frac_iv(x1).b])
        x1i = iv.mpf(frac_iv(x1).b)
        err = horner(MU_IV, x1i) * errB + NF * tailN(x1i)
        xm = frac_iv((x0 + x1) / 2)
        v1 = horner(Gc, xi)
        v2 = horner(Gc, xm) + horner(dGc, xi) * (xi - xm)
        low = max(lo(v1), lo(v2)) - hi(err)
        if low > 0:
            nx += 1
            vmin = low if vmin is None else min(vmin, low)
            lines.append(f"PASS box b in [{float(b0):.8f}, {float(b1):.8f}] beta^2 in [{float(x0):.8f}, {float(x1):.8f}]  FF >= {mp.nstr(low, 6)}")
        else:
            if depth >= 12 or nx + len(stack) > 4000:
                return False, [], None
            xmid = (x0 + x1) / 2
            stack.append((x0, xmid, depth + 1))
            stack.append((xmid, x1, depth + 1))
    return True, lines, vmin


def main():
    out = sys.stdout
    n0 = 8 * 64                               # initial b-cells
    cells = [(TAU2 * Fr(i, n0), TAU2 * Fr(i + 1, n0)) for i in range(n0)]
    nbox = 0
    ncell = 0
    vmin = None
    ok = True
    maxdepth = 0
    depth = 0
    with Pool(4) as pool:
        while cells and depth <= 12:
            res = pool.map(certify_cell, cells, chunksize=4)
            nxt = []
            for (b0, b1), (okc, lines, vm) in zip(cells, res):
                if okc:
                    ncell += 1
                    nbox += len(lines)
                    vmin = vm if vmin is None else min(vmin, vm)
                    out.write("\n".join(lines) + "\n")
                else:
                    bm = (b0 + b1) / 2
                    nxt += [(b0, bm), (bm, b1)]
            cells = nxt
            depth += 1
    if cells:
        ok = False
        for (b0, b1) in cells:
            out.write(f"FAIL b-cell [{float(b0):.8f}, {float(b1):.8f}]\n")
    out.write(f"SUMMARY C4: N = {N}, K = {K}; {ncell} b-cells, {nbox} boxes PASS; min certified lower bound of FF on R = {mp.nstr(vmin, 6)}\n")
    out.write(f"VERDICT C4 (FF > 0 on 0 <= b <= {float(TAU2)}, 0 <= beta <= b/{MK}): {'PASS' if ok else 'FAIL'}\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
