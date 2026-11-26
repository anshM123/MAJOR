"""rv_closed.py -- my implementation of the note's Sec.3 closed form ('uncut ideal angular action / w^3'), exact rationals.
   g[k,p](y): coefficient of z^k in (1-z)^(-p) P_{-1/2+i nu}((1+z)/(1-z)), y = nu^2, via the stated recurrence
   (validated against mpmath Taylor coefficients in rv6)."""
from fractions import Fraction as Fr
from functools import lru_cache
from rv_common import PSI2, PSI3, mom

@lru_cache(None)
def gpolys(p, N):
    """list of coefficient lists (ascending in y) of g[k,p], k=0..N"""
    out = [[Fr(1)]]
    prev = [Fr(0)]
    for k in range(N):
        cur = out[-1]
        nxt = [Fr(0)]*(len(cur)+1)
        for i, c in enumerate(cur):
            nxt[i] += (8*k*k+8*p*k+4*p-1)*c
            nxt[i+1] += -4*c
        for i, c in enumerate(prev):
            nxt[i] += -4*(k-1+p)**2*c
        nxt = [c/(4*(k+1)**2) for c in nxt]
        prev = cur
        out.append(nxt)
    return tuple(tuple(x) for x in out)

@lru_cache(None)
def moments2(N):
    return tuple(mom(PSI2, j) for j in range(2*N+3))

@lru_cache(None)
def moments3(N):
    return tuple(mom(PSI3, j) for j in range(2*N+3))

def AG(N, t):
    """A_k = 2 t^k int psi2 g[k,3](p^2/t) dp, G_k = t^k int psi3 g[k,4](p^2/t) dp, k = 0..N+1 (exact in t)."""
    g3 = gpolys(3, N+1); g4 = gpolys(4, N+1)
    m2 = moments2(N+1); m3 = moments3(N+1)
    A = [2*sum(c*t**(k-i)*m2[2*i] for i, c in enumerate(g3[k])) for k in range(N+2)]
    G = [sum(c*t**(k-i)*m3[2*i] for i, c in enumerate(g4[k])) for k in range(N+2)]
    return A, G

def closed_form(J, b, m2, lam2, AGcache=None):
    """sum_{k<=n} [u^k prod_{j<k} B_j / prod_{j<k}(1+B_j t u)] [A_k + u/(1+C t u)(t gamma A_k + (k+1)^2 A_{k+1}/2 + sigma G_k)]
    with t = 1/lam2, u = lam2/m2 (exact Fractions)."""
    t = Fr(1)/lam2; u = Fr(lam2)/m2
    n = J-b; C = J*(J+1); sig = 1 if b == 0 else -1
    A, G = AGcache if AGcache is not None else AG(n, t)
    tot = Fr(0); pref = Fr(1)
    fac = u/(1+C*t*u)
    for k in range(n+1):
        gam = Fr(-k*k+(8-4*b)*k+9-2*b*(b-1), 2)
        tot += pref*(A[k] + fac*(t*gam*A[k] + Fr((k+1)**2, 2)*A[k+1] + sig*G[k]))
        Bk = C-(b+k)*(b+k+1)
        pref = pref*u*Bk/(1+Bk*t*u)
    return tot
