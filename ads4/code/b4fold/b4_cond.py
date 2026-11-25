"""
b4_cond (BATCH4_FOLD.md §3.2, diagnostic): conditioning of hs_odd's bulk shape basis (hs_core._basis, column-scaled exactly as ShapeFit)
on the direct bulk nodes of the 2000 rows, per degree K = 4 ... 8, for the bulk's lower end XB = 0.025 (this note) and 0.05 (hs_odd),
and the twist-edge value of each degree (gk edge) -- BATCH5_FOLD §0 item 4 found the K = 6 basis ill-conditioned at J = 17, 19
(cond 6e7 / 1.7e7 against <= 4.6e5 where hs_odd was used before) and switched those states to the own-row route.
Usage: python3 b4_cond.py J,J,...   -> cache/b4_cond.txt (appended)
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b4_core as C
import numpy as np
import b4_hsodd as HO
H = C.H
lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


if __name__ == '__main__':
    for J in [int(v) for v in sys.argv[1].split(',')]:
        D = C.DG + J; rowm = H.HN.mphys(D, J); mtw = rowm/H.E.lam_of(C.DG); HO.XROW[0] = 1/mtw
        rK, rI, isdir = HO.rows(J, D)
        bad = H.refbad(J, D, HO.RTOL)
        for xb in (0.025, 0.05):
            HO.XB = xb
            xs = rI['xs']; sel = (xs >= xb - 1e-9) & (xs < HO.XCUT - 1e-9) & isdir
            out = []
            for K in (4, 5, 6, 7, 8):
                A = H._basis('+-', J, xs[sel], float(xs[sel].max()), K)
                sc = np.abs(A).max(0); sc[sc == 0] = 1
                cond = np.linalg.cond(A/sc)
                M = HO.Model(J, rK, rI, isdir, bad, bulk='dir', kb=K)
                lam = HO.value(J, mtw, M, 'gk')[0]
                out.append(f"K={K}: cond {cond:.1e}, tw {lam:+.4f}")
            pr(f"J={J} XB={xb}: {sel.sum()} bulk nodes | " + " | ".join(out))
    open(os.path.join(C.B4C, 'b4_cond.txt'), 'a').write("\n".join(lines) + "\n")
