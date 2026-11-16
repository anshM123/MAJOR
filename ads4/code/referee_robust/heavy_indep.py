"""Referee: independent floating-point scan of FF/N for M02/M01 (own quadrature + own Jacobi recurrences, coefficients typed
from ROBUST_FUNCTIONAL.md), plus exact values at small J (own hypergeometric sum)."""
import sys
import numpy as np
from fractions import Fraction as Fr
from low_indep import FUNS

def polymul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] += x * y
    return r
OMP3 = [1, -3, 3, -1]
def weights(name):
    x, y, a = FUNS[name]
    psi2 = polymul(OMP3, [0] + [Fr(c) for c in x][1:]) if False else polymul(OMP3, [Fr(0)] + [Fr(c) for c in x[1:]])
    psi3 = polymul(OMP3, [Fr(c) for c in y])
    return psi2, psi3, Fr(a)

def exact_parts(name, h, J, u):
    """exact B2, B3, FW at integer J, rational u via 2F1 series of D (P_J(1-2y)=2F1(-J,J+1;1;y), P^(0,8)_{J-4}(1-2y)=2F1(4-J,J+5;1;y))."""
    psi2, psi3, a = weights(name)
    sig = 1 if h == '++' else -1
    A_, B_ = (-J, J + 1) if h == '++' else (4 - J, J + 5)
    # D(1-2 u p^2) = sum_k c_k (u p^2)^k
    c = [Fr(1)]; k = 0
    while True:
        nxt = c[-1] * (A_ + k) * (B_ + k) / Fr((k + 1) ** 2)
        if nxt == 0: break
        c.append(nxt); k += 1
    mom = lambda w, j: sum(Fr(cc) / (i + j + 1) for i, cc in enumerate(w))
    I = lambda w: sum(ck * u ** kk * mom(w, 2 * kk) for kk, ck in enumerate(c))
    B2 = 2 * I(psi2) - u * I([Fr(0), Fr(0)] + psi2)
    B3 = sig * u * I(psi3)
    JJ = J * (J + 1); c4 = 1 if h == '++' else 41
    FW = a * (2 * JJ - c4 + sig) * u ** 3
    return B2, B3, FW

def scan(name, nq=700):
    psi2, psi3, a = weights(name)
    P2 = np.polynomial.polynomial.Polynomial([float(c) for c in psi2])
    P3 = np.polynomial.polynomial.Polynomial([float(c) for c in psi3])
    xg, wg = np.polynomial.legendre.leggauss(nq); pp = (xg + 1) / 2; ww = wg / 2
    us = np.unique(np.concatenate([np.logspace(-5, 0, 160), np.linspace(0.3, 1, 141)]))
    worst = {}
    for h in ('++', '+-'):
        sig = 1 if h == '++' else -1; c4 = 1 if h == '++' else 41
        best = (np.inf, None); cnt = 0; below = 0
        for u in us:
            Jm = int(min(2500, (200 / np.sqrt(u) - 1) / 2))
            x = 1 - 2 * u * pp * pp
            n_max = Jm if h == '++' else Jm - 4
            D = np.empty((n_max + 1, nq))
            if h == '++':
                D[0] = 1; D[1] = x
                for n in range(2, n_max + 1): D[n] = ((2 * n - 1) * x * D[n - 1] - (n - 1) * D[n - 2]) / n
                Js = np.arange(0, n_max + 1)
            else:
                al, be = 0.0, 8.0
                D[0] = 1; D[1] = (al + 1) + (al + be + 2) * (x - 1) / 2
                for n in range(2, n_max + 1):
                    c = 2 * n + al + be
                    D[n] = ((c - 1) * (c * (c - 2) * x + al * al - be * be) * D[n - 1] - 2 * (n + al - 1) * (n + be - 1) * c * D[n - 2]) / (2 * n * (n + al + be) * (c - 2))
                Js = np.arange(0, n_max + 1) + 4
            B2 = D @ (ww * P2(pp) * (2 - u * pp * pp))
            B3 = sig * u * (D @ (ww * P3(pp)))
            FW = float(a) * (2 * Js * (Js + 1) - c4 + sig) * u ** 3
            if h == '++':
                sel = Js % 2 == 0; Js, B2, B3, FW = Js[sel], B2[sel], B3[sel], FW[sel]
            F = B2 + B3 + FW; N = np.abs(B2) + np.abs(B3) + np.abs(FW)
            r = F / N; i = np.argmin(r); cnt += r.size; below += int(np.sum(r < 0.18))
            if r[i] < best[0]: best = (r[i], dict(m=1 / np.sqrt(u), J=int(Js[i]), bhat=(2 * Js[i] + 1) * np.sqrt(u)))
        worst[h] = best
        print('%s %s: %d states, min FF/N = %.5f at %s ; states with ratio < 0.18: %d' % (name, h, cnt, best[0], best[1], below), flush=True)
    return worst

if __name__ == '__main__':
    for name in sys.argv[1:] or ['M02']:
        # exact small-J values at m = M (u = 1) and margin threshold delta* = FF/N
        for h, Js in (('+-', range(4, 12)), ('++', range(0, 12, 2))):
            for J in Js:
                B2, B3, FW = exact_parts(name, h, J, Fr(1))
                F = B2 + B3 + FW; N = abs(B2) + abs(B3) + abs(FW)
                if J < 8: print('%s %s J=%d u=1: FF/N = %.6f (exact %s)' % (name, h, J, float(F / N), F / N if J == 4 else ''))
        scan(name)
