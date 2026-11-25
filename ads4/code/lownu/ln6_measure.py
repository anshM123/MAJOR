"""Check of the reviewer's Markov bounds for the spin-Gamma (Bernstein) mass distribution (exact rational arithmetic).
Model (reviewer): M^2 = A V, A ~ Gamma(alpha, scale 2), V ~ Gamma(gamma, scale 2) independent, after the M^{-2h} reweighting
alpha = lam - h, gamma = lam + J + 1/2 - r - h, lam = (tau-1)/2 (LARGESPIN_L L2: (-w^2)^{-(tau-1)/2}).
E[M^{-2n}] = prod_{k=1}^n 1/(4 (alpha-k)(gamma-k));  P(M < m0) <= min_n m0^{2n} E[M^{-2n}]  (Markov)."""
from fractions import Fraction as Fr
import math
from ln_common import guard
def bound(tau, J, r, h, m0, nmax=None):
    lam = Fr(tau - 1, 2); al = lam - h; ga = lam + J + Fr(1, 2) - r - h
    best = (Fr(1), 0); val = Fr(1)
    nmax = nmax or int(min(al, ga)) - 1
    for n in range(1, nmax):
        val = val * Fr(m0 * m0) / (4 * (al - n) * (ga - n))
        if val < best[0]: best = (val, n)
    return best
out = []
for tau in (1000, 2000):
    for m0 in (600, 700, 800, 900):
        for J in (0, 50, 489, 5000):
            for r in (0, 4):
                v, n = bound(tau, J, r, 7, m0)
                out.append('tau=%d m0=%d J=%d r=%d: P(M<m0) <= %.3e (n=%d)' % (tau, m0, J, r, float(v), n))
txt = '\n'.join(out); print(txt); open(guard('cache/ln6_measure.txt'), 'w').write(txt + '\n')
