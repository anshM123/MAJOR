"""
Same exact Mack polynomial P_{Delta,J}(s,t) as mack.mack_P_fast (Costa--Goncalves--Penedones, conventions of mack.py), built with FLINT's
fmpq_mpoly (C-level exact rational bivariate polynomials) instead of Python Fractions.  Identical algorithm and output layout
({(i, j): Fraction} = coefficient of s^i t^j), so skelexact3.Laplace3 can use it unchanged.  Checked against mack_P_fast (h13 self-test:
`python3 mack_flint.py`).  Delta should be a binary float with a short mantissa (e.g. a multiple of 1/8) to keep the rationals small.
"""
from fractions import Fraction as Fr
from math import factorial
import functools
import flint

CTX = flint.fmpq_mpoly_ctx.get(('s', 't'), 'lex')
S, T = CTX.gens()


def _q(x):
    x = Fr(x)
    return flint.fmpq(x.numerator, x.denominator)


def _rf(x, k):
    r = Fr(1)
    for q in range(k):
        r *= (x + q)
    return r


def _poch(lin, k):
    out = CTX.constant(1)
    for q in range(k):
        out = out*(lin + q)
    return out


@functools.lru_cache(None)
def mack_P_flint(Delta, J):
    D = Fr(Delta)
    h = Fr(3, 2)
    hp = D; hm = 2*h - D
    s1 = -S/2
    s2 = (S + T)/2
    a = 1 - (hp + J)/2; b = 1 - (hm + J)/2
    den = _rf(hp - 1, J)*_rf(hm - 1, J)
    pa = [_poch(s1, k) for k in range(J + 1)]
    pb = [_poch(s2, k) for k in range(J + 1)]
    fact = [factorial(i) for i in range(J + 1)]
    ra = [_rf(a, i) for i in range(J + 1)]; rb = [_rf(b, i) for i in range(J + 1)]
    pbb, paa = {}, {}
    P = CTX.constant(0)
    for r in range(J//2 + 1):
        aJr = Fr((-1)**r*factorial(J), 2**(2*r)*factorial(r)*factorial(J - 2*r))/_rf(h + J - 1 - r, r)
        pre = aJr*Fr(2)**(J + 2*r)*factorial(J - 2*r)/den
        l1 = _q((hp - J)/2) - T/2
        l2 = _q((hm - J)/2) - T/2
        pr = _poch(l1, r)*_poch(l2, r)
        K = J - 2*r; L = J - r
        inner = CTX.constant(0)
        for k13 in range(K + 1):
            for k24 in range(K + 1 - k13):
                rest = K - k13 - k24
                comb = CTX.constant(0)
                for k14 in range(rest + 1):
                    k23 = rest - k14
                    c = Fr((-1)**(k13 + k24), fact[k13]*fact[k14]*fact[k23]*fact[k24])
                    c *= ra[L - (k13 + k14)]*ra[L - (k23 + k24)]*rb[L - (k13 + k23)]*rb[L - (k14 + k24)]
                    if c == 0:
                        continue
                    key = (min(k14, k23), max(k14, k23))
                    if key not in pbb:
                        pbb[key] = pb[k14]*pb[k23]
                    comb = comb + pbb[key]*_q(c)
                if not comb.is_zero():
                    key = (min(k13, k24), max(k13, k24))
                    if key not in paa:
                        paa[key] = pa[k13]*pa[k24]
                    inner = inner + paa[key]*comb
        P = P + pr*inner*_q(pre)
    out = {}
    for mon, c in zip(P.monoms(), P.coeffs()):
        out[(int(mon[0]), int(mon[1]))] = Fr(int(c.p), int(c.q))
    return out


if __name__ == '__main__':
    import sys, time
    sys.path.insert(0, __file__.rsplit('/', 1)[0])
    import mack
    for J in (0, 1, 2, 5, 8, 13):
        for D in (900.25, 301.125):
            A = mack.mack_P_fast(D, J); B = mack_P_flint(D, J)
            ok = (set(A) == set(B)) and all(A[k] == B[k] for k in A)
            print('J', J, 'Delta', D, 'identical to mack_P_fast:', ok, flush=True)
    for J in (20, 30):
        t = time.time(); mack_P_flint(900.25, J); print('J', J, 'flint %.1fs' % (time.time() - t), flush=True)
