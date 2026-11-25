"""ANGULAR_GO Sec. 1: what the celestial kernels say about slots.
(a) j'=0 and j'=1 Fourier ratios k(m)/k(2) (kb_celestial Fourier sums, read-only import) = CDH Q_{m-2}, (d,l) = (1,0), m = 3..6;
    parity k^(0)(-m) = k^(0)(m), k^(1)(-m) = -k^(1)(m)  [DV].
(b) slot bookkeeping: heavy (J, lambda, h) -> (|m|, channel s); T family -> Einstein (m=0), WxW (|m|=4, s=-1)  [P].
(c) flat-regime alternation: sign Q_j(nu^2) = (-1)^j beyond the largest zero; largest zero vs j  [NE]; CDH weight ~ e^{-pi nu}  [P].
Usage: python3 -B ag1_kernels.py   (~1-2 min)  -> cache/ag1_kernels.txt"""
import sys, os
import mpmath as mp
from ag_common import Qrow, write, HERE
sys.path.insert(0, os.path.join(HERE, '..', 'kb2a4'))
from kb_celestial import Ssum  # read-only

out = []
def pr(s):
    print(s, flush=True); out.append(s)

mp.mp.dps = 30
pr('# ag1_kernels')
pr('(a) Fourier ratios k(m)/k(2) vs Q_{m-2}(nu^2), (d,l)=(1,0); N = 2000 modes')
worst = 0
for nu in ('0.7', '1.9', '3.3'):
    x = mp.mpf(nu)**2
    Q = Qrow(1, 0, x, 4)
    for jp in (0, 1):
        k2 = Ssum(nu, 2, jp, 2000)
        for m in (3, 4, 5, 6):
            r = Ssum(nu, m, jp, 2000)/k2
            dev = abs(r/Q[m - 2] - 1); worst = max(worst, dev)
            rm = Ssum(nu, -m, jp, 2000)/k2
            par = abs(rm/r - (1 if jp == 0 else -1))
            pr('   nu=%s jp=%d m=%d: k(m)/k(2)=%s  Q=%s  rel dev %.1e | k(-m)/k(m) parity dev %.1e' %
               (nu, jp, m, mp.nstr(r.real, 12), mp.nstr(Q[m - 2], 12), float(dev), float(par)))
pr('   max rel deviation %.1e (truncation-limited)' % float(worst))

pr('(b) slots.  heavy family spin J, helicity lambda in [-J,J], external h = +-2; harmonic m = lambda - h, |m| >= 2 (Ward);')
pr('    k=3 dual: d_-- = -d_++, kernel odd  => slot value  A2(|m|) + s A3(|m|),  s = sgn(h m)')
for M in range(2, 9):
    pops = {+1: [], -1: []}
    for h in (2, -2):
        for lam in (h + M, h - M):
            s = 1 if h*(lam - h) > 0 else -1
            pops[s].append('J>=%d (lam=%+d,h=%+d)' % (abs(lam), lam, h))
    pr('   |m|=%d: s=+1 from %s ; s=-1 from %s' % (M, pops[1][0], min(pops[-1], key=lambda t: int(t[3:t.index(' ')]))))
pr('   T family (lambda = +-2 only): Einstein lambda=h -> m=0 (no heavy slot); WxW lambda=-h -> |m|=4, s=-1 (heavy slot (4,-): J>=2)')

pr('(c) largest zero of Q_j(nu^2) and the sign beyond it')
for d, l in ((1, 0), (0, 0)):
    for j in (2, 4, 8, 16, 32, 64):
        # bisection on a grid in nu
        grid = [mp.mpf(i)/20 for i in range(1, 20*(3*j + 20))]
        vals = [Qrow(d, l, g**2, j)[j] for g in grid]
        last = max(i for i in range(len(grid) - 1) if vals[i]*vals[i + 1] < 0)
        nu_z = grid[last]
        sgn = 1 if vals[-1] > 0 else -1
        pr('   (d,l)=(%d,%d) j=%2d: largest zero nu ~ %.2f (nu/j = %.3f); sign for nu above it = %+d = (-1)^j: %s' %
           (d, l, j, float(nu_z), float(nu_z/j), sgn, sgn == (-1)**j))
a, b, c = mp.mpf(3)/2, mp.mpf(1)/2, mp.mpf(1)/2
for nu in (5, 10, 20):
    w = abs(mp.gamma(a + 1j*nu)*mp.gamma(b + 1j*nu)*mp.gamma(c + 1j*nu)/mp.gamma(2j*nu))**2
    pr('   CDH weight (d,l)=(0,0) at nu=%d: w/(nu^2/cosh(pi nu)) = %s' % (nu, mp.nstr(w/(nu**2/mp.cosh(mp.pi*nu)), 10)))
write('ag1_kernels.txt', out)
