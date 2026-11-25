"""Shared helpers for ads4/code/angular2 (ANGULAR_GO.md).  Writes only inside ./cache (guard)."""
import os
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')


def guard(path):
    p = os.path.abspath(path)
    assert p.startswith(CACHE + os.sep), 'refusing to write outside angular2/cache: ' + p
    return p


def Qrow(d, l, x, jmax):
    """normalised continuous dual Hahn Q_j(x), j = 0..jmax, (a,b,c) = (2d+3/2, 1/2, l+1/2), x = nu^2, via the SP1 recurrence
    (k_{d} = 0, k_{d+1} = 1; Q_j = k_{j+d+1}).  mpmath, any precision."""
    k = {d: mp.mpf(0), d + 1: mp.mpf(1)}
    for m in range(d + 1, d + 1 + jmax):
        k[m + 1] = -((x + mp.mpf(1)/4 + 2*d*(d + 1) - 2*m*(m + l))*k[m] + (m - 1 - d)*(m + l - 1 - d)*k[m - 1]) / ((m + 1 + d)*(m + l + 1 + d))
    return [k[j + d + 1] for j in range(jmax + 1)]


def write(name, lines):
    open(guard(os.path.join(CACHE, name)), 'w').write('\n'.join(lines) + '\n')
