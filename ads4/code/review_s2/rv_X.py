"""rv_X.py [DV]: Lemma X / helicity (rt_pp.PP) at the right edge nu = Lam + iy of the contour, which carries 0.32 F of the [3.65, 3.75]
certificate.  The right edge uses PP(sigma_top, Lam - 1, Lam_up, d0) for sigma' in [0, sigma_top]; we compare with exact
max |S^(sigma' + i tau)| over tau in {Lam - 1, Lam - 0.5, Lam_up}, sigma' in a grid of [0, sigma_top], integer a in [998, 1200], both b,
n at the cell edges L = 3.65, 3.75 (fastS)."""
import sys, math
sys.path.insert(0, '/home/user/MAJOR/ads4/code/review_s2/ws/code/analytic2')
import common2
import rt_pp as PPm
import rt_cert as RC
from fastS import S_vec
La, Lb = (float(sys.argv[1]), float(sys.argv[2])) if len(sys.argv) > 2 else (3.65, 3.75)
d_hi = math.exp(-La); d0 = d_hi*1.0001
LAMf, LAMup = RC.LAMf, RC.LAMup


def dval(n, a, s):
    N = n + a
    return a*(2*n + a + s)/(N*(N + s))


def n_at(a, s, dt):
    lo, hi = 1, 10**12
    while hi - lo > 1:
        m = (lo + hi)//2
        if dval(m, a, s) > dt: lo = m
        else: hi = m
    return hi


worst = 1e9
for stop in (0.5, 1.5, 2.5, 3.5):
    bnd = PPm.PP(stop, LAMf - 1, LAMup, d0)
    taus = [LAMf - 1, LAMf - 0.5, LAMup]
    sigs = [stop*k/4 for k in range(5)]
    exm = -1e9; at = None
    for a in (998, 999, 1000, 1002, 1005, 1010, 1020, 1050, 1100, 1200):
        for b in (0, 4):
            s = 2*b + 1
            for L in (La, Lb):
                n = n_at(a, s, math.exp(-L))
                dd = dval(n, a, s)
                als = [complex(sg, t) for sg in sigs for t in taus]
                S, lg = S_vec(n, a, s, als)
                for k, al in enumerate(als):
                    ex = math.log(abs(S[k])) + lg[k] + al.real*math.log(dd)
                    if ex > exm: exm = ex; at = (a, b, n, al)
    gap = math.log(bnd) - exm
    worst = min(worst, gap)
    print(f'PP(sigma<={stop}, tau in [Lam-1, Lam]) : log bound {math.log(bnd):.3f}, max exact {exm:.3f} at {at}, gap {gap:.3f}', flush=True)
print('min gap', worst)
