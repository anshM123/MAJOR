"""
lsw/lsw_r1rank.py (LARGESPIN_W.md §9): test of the R1 sufficient condition (§3.4): does the leading-r EE kernel factorise as
    L_nu[w1, w3, z2, z4](B) = T(pols; B) * l(nu, B)      (one polarisation tensor)?
L_nu = r^0 Cauchy coefficient of K^EE/(N r^2) (EESTEP §1.4), computed with eestep's own code (ee6_lrprofile.L, read-only import).
Test: at fixed B, the matrix M[pol set, nu] has rank 1 iff the factorisation holds.  Singular values printed.
Usage: python3 lsw_r1rank.py -> cache/r1rank.txt
"""
import os, sys
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'eestep'))
import numpy as np
import ee6_lrprofile as E6          # read-only

out = []
pr = lambda s: (print(s, flush=True), out.append(s))
rng = np.random.default_rng(11)
PS = [tuple(rng.uniform(0, 2*np.pi, 4)) for _ in range(7)]
NUS = [0.5, 1.5, 3.0, 6.0, 10.0]
for B in (0.3, 1.0, 3.0):
    M = np.array([[E6.L(nu, B, pa)[0] for nu in NUS] for pa in PS])
    sv = np.linalg.svd(M, compute_uv=False)
    P = np.array([E6.L(nu, B, PS[0])[1] for nu in NUS])
    pr('B=%.2f  singular values / s0: %s' % (B, ' '.join('%.2e' % (x/sv[0]) for x in sv)))
    pr('        row0 L/P_l: %s' % ' '.join('%+.4e%+.4ej' % (z.real, z.imag) for z in M[0]/P))
open(os.path.join(HERE, 'cache', 'r1rank.txt'), 'w').write('\n'.join(out) + '\n')
