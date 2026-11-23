r"""
rh_common (R4ARC_HEFT.md): shared helpers.  Read-only imports of LIGHTSIDE_DIST's machinery (resultq/lsd_core, lsd_mellin, lsd_pos,
og_free); nothing there is modified.

Object: the swapped-form functional Psi_k at FIXED celestial angles (a, b) on the product correlator T' = sigma x T, split by the Wick
class c of the T-factor (og_free: connected free classes 1234, 1243, 1324; disconnected/MFT classes s, t, u).  For a sigma-part with
Mellin amplitude M(d12, d13) (measure Gamma(d12)^2 Gamma(d13)^2 Gamma(d14)^2, d14 = 2 - d12 - d13, Delta_sigma = 2), the term-wise
(s-channel OPE) value is
    Omega[M x c](a, b) = (1/2pi) int dy13 Gamma(d13)^2 sum_{n>=0} Res_{d12=-n}[Gamma(d12)^2 Gamma(d14)^2 M K_c(d12, d13)]
(+ single-trace residues of M if M has s-channel poles), d13 = c2 + i y13, exactly as lsd_core.residues_d12, but for an ARBITRARY
amplitude given as a callable returning (M, dM/dd12) at fixed d13.
"""
import os, sys
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
RQ = os.path.join(HERE, '..', 'resultq')
sys.path.insert(0, RQ)
import numpy as np
from scipy.special import loggamma, digamma

import lsd_core as C      # noqa: E402  (read-only)
import lsd_pos as P       # noqa: E402
import og_free as F       # noqa: E402

CACHE = os.path.join(HERE, 'cache')
GRID = dict(s0=1e-9, m=16, q=1.8)        # LIGHTSIDE_DIST / F_NORM production grid


def setup(k, phi1, phi3, p2, p4, classes=C.CLS, grid=GRID):
    return C.Setup(k, phi1, phi3, p2, p4, classes=classes, **grid)


def residues(S, Mfun, d13, N):
    """left double-pole residues at d12 = -n, n = 0..N, of Gamma(d12)^2 Gamma(d14)^2 M K_c at fixed d13 (no Gamma(d13)^2).
    Mfun(d12_array, d13) -> (M, dM/dd12).  Returns {c: (N+1,) array}."""
    nf = np.arange(N + 1).astype(float)
    d12 = -nf.astype(complex)
    d14 = 2 - d12 - d13
    Kd = S.K(d12, d13, deriv=True)
    g = np.exp(2*loggamma(d14) - 2*loggamma(nf + 1))
    psi14 = digamma(d14); psin = digamma(nf + 1)
    M, dM = Mfun(d12, d13)
    out = {}
    for c, (K, dK) in Kd.items():
        out[c] = g*((-2*psi14*M + dM)*K + M*dK + 2*psin*M*K)
    return out


def d13_line(c2=2/3, Y=8.0, h=0.2):
    y = np.arange(-Y, Y + h/2, h)
    return c2 + 1j*y, h


def poly_D(J):
    """top Regge monomial D^J/4, D = d12 - d14 = 2 d12 + d13 - 2 (r4 of LIGHTSIDE_DIST for J = 4)."""
    def f(d12, d13):
        D = 2*d12 + d13 - 2
        return D**J/4, 2*J*D**(J - 1)/4
    return f


def poly_q(j):
    """q_j = d12^j + d13^j + d14^j (crossing-symmetric EFT monomial sum)."""
    def f(d12, d13):
        d14 = 2 - d12 - d13
        if j == 0:
            return 3*np.ones_like(d12), np.zeros_like(d12)
        return d12**j + d13**j + d14**j, j*d12**(j - 1) - j*d14**(j - 1)
    return f


def tail(r, n0):
    """power-law tail estimate of sum_{n > N} r_n from the residues on [n0, N] (LIGHTSIDE_DIST / fn_spin convention)."""
    n = np.arange(n0, len(r)); yv = r[n0:]
    if np.all(np.abs(yv) < 1e-200):
        return 0.0, 0.0
    pp = -np.polyfit(np.log(n), np.log(np.abs(yv) + 1e-300), 1)[0]
    A = np.sum(yv*n**pp)/len(n)
    if pp <= 1.05:
        return pp, np.nan
    return pp, A*(len(r) - 0.5)**(1 - pp)/(pp - 1)


def values(Sab, Sba, k, amps, N=160, line=None, nt=80):
    """symmetrised term-wise values per class: V = E(a,b) + (-1)^k E(b,a), E = (1/2pi) int dy13 Gamma(d13)^2 sum_n Res.
    amps: {name: Mfun}.  Returns {(name, c): (value_without_tail, tail, p)} and the raw residue arrays."""
    d13s, h = line if line is not None else d13_line()
    acc = {}
    for d13 in d13s:
        w = np.exp(2*loggamma(d13))*h/(2*np.pi)
        for o, S in (('ab', Sab), ('ba', Sba)):
            sg = 1.0 if o == 'ab' else (-1.0)**k
            for name, f in amps.items():
                r = residues(S, f, d13, N)
                for c, v in r.items():
                    acc[(name, c)] = acc.get((name, c), 0) + sg*w*v
    out = {}
    for key, r in acc.items():
        pp, tl = tail(r, nt)
        out[key] = (np.sum(r), tl, pp)
    return out, acc
