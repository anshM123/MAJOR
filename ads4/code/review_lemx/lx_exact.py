"""lx_exact.py [DV]: Lemma X (rt_pp.PP0 = min_c cross(sigma, Lam-1, Lam_up, c), uniform a >= 998) against exact b = 0 multipliers
|S^(sigma + i tau)| on the right edge of the [3.65, 3.75] contour (tau = Lam; the bound covers tau in [Lam - 1, Lam_up]).
Stage 1: double-precision scan (lx_common.S_hat_vec) over a in [998, 4500] (incl. non-integer a near 998, where y_max is largest),
         n at L in {3.65, 3.70, 3.75}, sigma in [0, 5.5], tau in {Lam-1, Lam-0.5, Lam}; also an all-spin sweep at a = 998.
Stage 2: Arb confirmation (12000 and 16000 bits) of the 8 tightest points; b = 4 (s = 9, fastS) vs PP at a few states."""
import math, sys, time
sys.path.insert(0, '/home/user/MAJOR/ads4/code/analytic'); sys.path.insert(0, '/home/user/MAJOR/ads4/code/analytic2')
import rt_pp
from lx_common import S_hat_vec, S_hat_arb, n_at, d_of, LAMf
SIG = [0.0, 0.25, 0.5, 1.0, 1.5, 2.5, 3.5, 4.5, 5.5]
TAU = [LAMf - 1, LAMf - 0.5, LAMf]
PB = {s: math.log(rt_pp.PP0(s, LAMf - 1, LAMf)) for s in SIG}
print('log PP0(sigma, [Lam-1, Lam]) =', {s: round(v, 3) for s, v in PB.items()}, flush=True)
AS = [998.0, 998.1, 998.25, 998.4, 998.49, 998.5, 998.75, 999.0, 999.5, 1000.0, 1002.0, 1005.0, 1010.0, 1020.0, 1050.0, 1100.0,
      1200.0, 1500.0, 2000.0, 3000.0, 4500.0]
rows = []
t0 = time.time()
als = [complex(s, t) for s in SIG for t in TAU]
for a in AS:
    for L in (3.65, 3.70, 3.75):
        n = n_at(a, math.exp(-L))
        v = S_hat_vec(n, a, als)
        for i, al in enumerate(als):
            rows.append((PB[al.real] - v[i], a, n, L, al.real, al.imag, v[i]))
    best = min(r for r in rows if r[1] == a)
    print(f'a={a}: min gap {best[0]:.3f} at n={best[2]} L={best[3]} sigma={best[4]} tau={best[5]:.3f} (log|S^|={best[6]:.3f}) [{time.time()-t0:.0f}s]', flush=True)
rows.sort()
print('\nper sigma: min gap over the cell scan')
for s in SIG:
    r = min(x for x in rows if x[4] == s)
    print(f'  sigma={s}: log bound {PB[s]:.3f}, max log|S^| {r[6]:.3f} (a={r[1]}, n={r[2]}), gap {r[0]:.3f}')
print('\nall-spin sweep at a = 998, tau = Lam (Lemma X claims uniformity in n):')
for n in (1001, 5000, 20000, 80000, 300000, 1000000):
    v = S_hat_vec(n, 998.0, [complex(s, LAMf) for s in SIG])
    print(f'  n={n} (L={-math.log(d_of(n, 998.0)):.2f}): min gap {min(PB[s] - v[i] for i, s in enumerate(SIG)):.3f}; log|S^| = {[round(float(x), 3) for x in v]}', flush=True)
print('\nStage 2: Arb confirmation of the tightest points')
seen = set()
for r in rows:
    key = (r[1], r[2], r[4], r[5])
    if key in seen: continue
    seen.add(key)
    if len(seen) > 8: break
    lo, hi, rel = S_hat_arb(r[2], r[1], r[4], r[5], 12000)
    lo2, hi2, rel2 = S_hat_arb(r[2], r[1], r[4], r[5], 16000) if len(seen) <= 2 else (lo, hi, rel)
    print(f'  a={r[1]} n={r[2]} sigma={r[4]} tau={r[5]:.4f}: double {r[6]:.6f}, Arb12000 [{lo:.6f}, {hi:.6f}], Arb16000 [{lo2:.6f}, {hi2:.6f}];'
          f' bound {PB[r[4]]:.4f}; certified gap >= {PB[r[4]] - hi:.4f}', flush=True)
print('\nb = 4 (s = 9) vs PP (helicity on top of Lemma X), a few states')
from fastS import S_vec
for a in (998.0, 998.49, 1000.0, 1100.0):
    for L in (3.65, 3.75):
        N_ = None
        dt = math.exp(-L)
        lo_, hi_ = 1, 10**9
        dv = lambda n: a*(2*n + a + 9)/((n + a)*(n + a + 9))
        while hi_ - lo_ > 1:
            m = (lo_ + hi_)//2
            if dv(m) > dt: lo_ = m
            else: hi_ = m
        n = hi_; dd = dv(n)
        als4 = [complex(s, LAMf) for s in SIG]
        S, lg = S_vec(n, a, 9, als4)
        ex = [math.log(abs(S[i])) + lg[i] + SIG[i]*math.log(dd) for i in range(len(SIG))]
        gaps = [math.log(rt_pp.PP(s, LAMf - 1, LAMf, math.exp(-3.65)*1.0001)) - ex[i] for i, s in enumerate(SIG)]
        print(f'  a={a} L={L} n={n}: min gap (PP vs exact b=4) {min(gaps):.3f}', flush=True)
