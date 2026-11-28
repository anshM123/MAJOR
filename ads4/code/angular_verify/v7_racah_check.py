"""v7: independent exact-rational check of the claimed closed-form (Racah / 4F3) dual on a finite state grid.
Dual functional (rule multipliers mu_k, k = 0..n):
    mu_k binom(k+3,3)^2 = A (-n)_k (n+s)_k / (k! (s+1)_k),   A = -W_4(g+4)
so that f_b(J,D) = A * 4F3-type sum  sum_k (-n)_k (n+s)_k/(k!(s+1)_k) prod_{j<k} (C-y_{b+j})/(T-y_{b+j}).
Feasibility: W_4 + f_4 >= 0 on every b=4 state, W_0 - f_0 >= 0 on every b=0 state (heavy: twist >= g).
Bound: |c|^2 <= 3 / L,  L = -mu_2 = -A (-n)_2 (n+s)_2 / (2 (s+1)_2) / 100.
Also reports the margin on LIGHT twist-2 states (A2' coefficient).   usage: python3 v7_racah_check.py g n [s]"""
import sys
from fractions import Fraction as Fr
from math import sqrt
from common import W, allowed

g = int(sys.argv[1]); n = int(sys.argv[2]); s = int(sys.argv[3]) if len(sys.argv) > 3 else 0
A = -W(4, Fr(g + 4))
coef = [Fr(1)]
for k in range(n):
    coef.append(coef[-1] * (-n + k) * (n + s + k) / ((k + 1) * (s + 1 + k)))


def F(b, J, D):
    C = J * (J + 1); T = (D - 1) * (D - 2)
    tot = Fr(0); p = Fr(1)
    for k in range(n + 1):
        tot += coef[k] * p
        if k + 1 > J - b:
            break
        y = (b + k) * (b + k + 1)
        p = p * (C - y) / (T - y)
    return A * tot


def margin(b, J, D):
    f = F(b, J, D)
    w = W(b, D)
    return (w + f) / w if b == 4 else (w - f) / w


L = -A * coef[2] / 100
print(f"g={g} n={n} s={s}: claimed bound tau g^3 = {4*sqrt(3/float(L))*g**3:.2f}")
ts = [Fr(0), Fr(1, 4), Fr(1, 2), Fr(1), Fr(2), Fr(4), Fr(8), Fr(16), Fr(32), Fr(64), Fr(128), Fr(512), Fr(4096), Fr(10 ** 5)]
Js = list(range(0, 6 * g)) + [8 * g, 12 * g, 20 * g, 50 * g, 200 * g]
worst = (10.0, None)
for J in Js:
    for b in (0, 4):
        if not allowed(b, J):
            continue
        for t in ts:
            m = float(margin(b, J, Fr(J + g) + t))
            if m < worst[0]:
                worst = (m, (b, J, float(t)))
print(f"   heavy grid: min relative margin = {worst[0]:.4e} at (b,J,twist-g) = {worst[1]}  (>= 0 means feasible on grid)")
wl = (0.0, None)
for J in range(0, 4 * n):
    for b in (0, 4):
        if not allowed(b, J):
            continue
        D = Fr(J + 2)
        f = F(b, J, D)
        ex = (f if b == 4 else -f) - (Fr(3) / (1 / L * 1)) * 0   # excess of sigma f over 0
        if float(-ex) < wl[0]:
            pass
        val = float((f if b == 4 else -f))
        if val > wl[0]:
            wl = (val, (b, J))
print(f"   light twist-2 states: max sigma*f per unit |a_0|^2 = {wl[0]:.3e} at {wl[1]}  (units: endpoint budget 3; "
      f"compare L = {float(L):.3e})")
