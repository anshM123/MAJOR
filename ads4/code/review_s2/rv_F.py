"""rv_F.py [DV, adversarial]: the certificate's base bound PP_box (Lemma F / E' / Lemma X + helicity) on the actual certificate boxes
(unit a-box, tau-cell [T0, T0+1.1], L-cell) against exact |S^_n(sigma + i tau)| at integer corner/edge states of the box (both channels).
Exact values: dual-Hahn recurrence in n (fastS, complex128 with rescaling); the tightest cases are re-done in Arb (common2.S_rec).
Usage: python3 rv_F.py La Lb  (prints min gap log(bound) - log(exact) per box)."""
import sys, os, math
sys.path.insert(0, '/home/user/MAJOR/ads4/code/review_s2/ws/code/analytic2')
import common2
import rt2_kernel as RK
from fastS import S_vec
from flint import arb, acb, ctx

La, Lb = float(sys.argv[1]), float(sys.argv[2])
d_hi, d_lo = math.exp(-La), math.exp(-Lb)
T0s = [float(x) for x in (sys.argv[3].split(',') if len(sys.argv) > 3 else ['1003', '1020', '1031', '1040'])]
alos = [float(x) for x in (sys.argv[4].split(',') if len(sys.argv) > 4 else ['998', '1005', '1020', '1030'])]


def dval(n, a, s):
    N = n + a
    return a*(2*n + a + s)/(N*(N + s))


def n_range(a, s):
    """smallest/largest integer n >= 1001 - b with d in [d_lo, d_hi] (d decreasing in n)."""
    lo, hi = 1, 10**9
    while hi - lo > 1:
        m = (lo + hi)//2
        if dval(m, a, s) > d_hi: lo = m
        else: hi = m
    nmin = hi
    lo, hi = 1, 10**9
    while hi - lo > 1:
        m = (lo + hi)//2
        if dval(m, a, s) >= d_lo: lo = m
        else: hi = m
    nmax = lo
    return nmin, nmax


worst = 1e9
for T0 in T0s:
    T1 = T0 + 1.1
    for alo in alos:
        ahi = alo + 1
        for sig in (2.5, 3.5):
            bnd = RK.PP_box(sig, T0, T1, alo, ahi, d_lo, d_hi)
            b0 = RK.PP0_box(sig, T0, T1, alo, ahi, d_lo, d_hi)
            taus = [T0, T0 + 0.55, T1]
            ex_max = -1e9; arg = None
            for a in (int(alo), int(ahi)):
                for b in (0, 4):
                    s = 2*b + 1
                    nmin, nmax = n_range(a, s)
                    for n in (nmin, nmin + 1, (nmin + nmax)//2, nmax - 1, nmax):
                        S, lg = S_vec(n, a, s, [complex(sig, t) for t in taus])
                        dd = dval(n, a, s)
                        for k, t in enumerate(taus):
                            ex = math.log(abs(S[k])) + lg[k] + sig*math.log(dd)
                            if ex > ex_max:
                                ex_max = ex; arg = (a, b, n, t)
            gap = math.log(bnd) - ex_max
            worst = min(worst, gap)
            print(f'L[{La},{Lb}] tau[{T0},{T1}] a[{alo},{ahi}] sig={sig}: log PP_box={math.log(bnd):.3f} (b=0 part {math.log(b0):.3f})'
                  f'  max exact={ex_max:.3f} at (a,b,n,tau)={arg}  gap={gap:.3f}', flush=True)
print('min gap', worst)
