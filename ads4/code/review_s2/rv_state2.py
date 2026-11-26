"""rv_state2.py [exact]: state relations on the cell L in [3.65, 3.75] (and the whole range L >= 3.65), exact rationals, both channels.
Proof (exact algebra, see REVIEW_SLIVER2 §5): d = a(2N - a + s)/(N(N+s)), x = a/N:
  a/(Nd) = (N+s)/(2N-a+s) >= 1/2;  x(2-x) <= d(1+s/N) => x <= 0.01314 (d <= e^-3.65, s/N <= 9/1995) => N >= 75951, s/N <= 1.185e-4
  => x <= d(1+s/N)/(2-x) <= 0.50337 d;  D d^2/a^2 = (2-x+s/N)^2/(1+s/N) >= (2-0.01314)^2/1.0001185 = 3.9471;  N/a >= 76.
The scan below checks the extremes at the cell edges for a in {998..10^10} + random, minimal and maximal n in the cell."""
import math, random
from fractions import Fraction as Fr
random.seed(5)
Lmin, Lmax = 3.65, 3.75
d_hi = math.exp(-Lmin)


def dexact(n, a, s):
    N = n + a
    return Fr(a*(2*N - a + s), N*(N + s))


def n_for(a, s, dtarget):
    lo, hi = 1, 10**14
    while hi - lo > 1:
        m = (lo + hi)//2
        if float(dexact(m, a, s)) > dtarget: lo = m
        else: hi = m
    return hi          # smallest n with d <= dtarget


# analytic constants
x1 = 1 - math.sqrt(1 - d_hi*(1 + 9/1995))
Nmin = 998/x1; e = 9/Nmin
print(f'proof constants at d = e^-3.65 = {d_hi:.6f}: x <= {x1:.6f}, N >= {Nmin:.0f}, s/N <= {e:.4e}, '
      f'a/(Nd) <= {(1 + e)/(2 - x1):.5f}, Dd^2/a^2 >= {(2 - x1)**2/(1 + e):.5f}, N/a >= {1/x1:.1f}')
worst = dict(rmax=0, rmin=9, Dmin=99, Nmin=1e99)
for b in (0, 4):
    s = 2*b + 1
    for a in [998, 999, 1000, 1031, 1500, 3000, 10**4, 10**6, 10**8, 10**10] + [random.randint(998, 10**7) for _ in range(60)]:
        for L in (Lmin, Lmin + 0.03, Lmin + 0.07, Lmax):
            n = n_for(a, s, math.exp(-L))
            for nn in (n, n + 1, n - 1):
                if nn < 1001 - b: continue
                d = dexact(nn, a, s)
                if not (math.exp(-Lmax) <= float(d) <= d_hi): continue
                N = nn + a; D = N*(N + s)
                r = Fr(a, N)/d; Dd = D*d*d/(a*a)
                worst['rmax'] = max(worst['rmax'], float(r)); worst['rmin'] = min(worst['rmin'], float(r))
                worst['Dmin'] = min(worst['Dmin'], float(Dd)); worst['Nmin'] = min(worst['Nmin'], N/a)
print('scan L in [3.65, 3.75]:', {k: round(v, 6) for k, v in worst.items()})
assert worst['rmin'] >= 0.49 and worst['rmax'] <= 0.51 and worst['Dmin'] >= 3.9 and worst['Nmin'] >= 40
print('relations 0.49d <= a/N <= 0.51d, Dd^2 >= 3.9a^2, N >= 40a: hold on the scan (and by the proof above)')
