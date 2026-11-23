"""Target 1a: the magic-block de-shift identity for every b-shifted (0, b) block (PW_IMAGES (2.1)),
       2 eta (1+r^2) Ghat^1_b(nu) = (1 + i/2nu) Ghat^0_b(nu+i) + (1 - i/2nu) Ghat^0_b(nu-i),
checked with the independent recursion of blocks2.py.  Output: cache/t1_magic.txt
"""
import sys, time, random
import sympy as sp
import mpmath as mp
from blocks2 import magic_residual, recursion

out = open('cache/t1_magic.txt', 'w')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()

t0 = time.time()
# ---------------------------------------------------------------- (a) exact, l AND b symbolic
l, b = sp.symbols('l b')
conv = lambda e: e
simp = lambda e: sp.cancel(sp.together(e))
for Nmax in (4, 5, 6):
    res, lhs, _ = magic_residual(l, b, Nmax, conv, simp)
    nz = sum(1 for d in res for c in d.values() if c != 0)
    ncoef = sum(1 for d in lhs for c in d.values() if simp(c) != 0)
    log(f'(a) exact, l and b symbolic, N <= {Nmax}: nonzero residual coefficients = {nz} (of {ncoef} nonzero LHS coefficients)  [{time.time()-t0:.0f}s]')

# exact, l symbolic, integer b = -4..4, N <= 8
for bb in range(-4, 5):
    res, lhs, _ = magic_residual(l, sp.Integer(bb), 8, conv, simp)
    nz = sum(1 for d in res for c in d.values() if c != 0)
    ncoef = sum(1 for d in lhs for c in d.values() if simp(c) != 0)
    log(f'(a2) exact, l symbolic, b = {bb:+d}, N <= 8: nonzero residuals = {nz} of {ncoef} coefficients  [{time.time()-t0:.0f}s]')

# ---------------------------------------------------------------- (b) 50 digits, complex nu, N <= 40
mp.mp.dps = 50
convm = lambda e: mp.mpmathify(complex(e)) if e.is_number else None
def conv_mp(e):
    return mp.mpc(sp.re(e).evalf(60), sp.im(e).evalf(60)) if not e.is_real else mp.mpf(e.evalf(60))
random.seed(7)
nus = [mp.mpf('0.7'), mp.mpc('2.3', '0.2'), mp.mpc(random.uniform(0.2, 6), random.uniform(-0.4, 0.4)), mp.mpf(5)]
blist = [-4, -3, -2, -1, 0, 1, 2, 3, 4, mp.mpf('0.37'), 7, 10]
for nu in nus:
    lv = mp.mpf(-0.5) + 1j * nu
    row = []
    for bb in blist:
        res, lhs, _ = magic_residual(lv, bb, 40, conv_mp)
        scale = max(abs(c) for d in lhs for c in d.values())
        worst = max(abs(c) / max(1e-300, max(abs(x) for x in lhs[N].values()) if lhs[N] else 1)
                    for N, d in enumerate(res) for c in d.values())
        row.append((bb, worst))
    log(f'(b) 50 digits, nu = {mp.nstr(nu, 6)}, N <= 40: max |residual|/max|LHS at same N| per b: ' +
        ', '.join(f'b={mp.nstr(bb,3)}: {mp.nstr(w, 2)}' for bb, w in row) + f'  [{time.time()-t0:.0f}s]')

# ---------------------------------------------------------------- (c) controls: the identity is not trivially satisfied
nu = mp.mpc('2.3', '0.2'); lv = mp.mpf(-0.5) + 1j * nu
def worst_rel(res, lhs):
    return max(abs(c) / max(abs(x) for x in lhs[N].values()) for N, d in enumerate(res) for c in d.values() if lhs[N])
for bb in (0, 2, -3):
    r1, L1, _ = magic_residual(lv, bb, 12, conv_mp, cplus=2 * lv / (2 * lv + 1), cminus=2 * (lv + 1) / (2 * lv + 1))
    r2, L2, _ = magic_residual(lv, bb, 12, conv_mp, Dr=1)
    r3, L3, _ = magic_residual(lv, bb, 12, conv_mp, Dl=2, Dr=1)
    # b-dependent coefficients test: use the b = 0 blocks on the RHS for b != 0
    log(f'(c) controls b={bb}: swapped coefficients {mp.nstr(worst_rel(r1, L1), 3)}; RHS Delta=1 {mp.nstr(worst_rel(r2, L2), 3)}; '
        f'(Delta_L, Delta_R) = (2, 1) {mp.nstr(worst_rel(r3, L3), 3)}')
# mixing b between sides must fail (b-dependence of the blocks is real)
from blocks2 import eta_mult_basis
B1 = recursion(lv, 2, 1, 10, conv_mp); Bm = recursion(lv - 1, 0, 0, 10, conv_mp); Bp = recursion(lv + 1, 0, 0, 10, conv_mp)
diff = max(abs(B1[N].get(m, 0) - recursion(lv, 0, 1, 10, conv_mp)[N].get(m, 0)) for N in range(10) for m in B1[N])
log(f'(c2) the blocks do depend on b: max |B^(b=2) - B^(b=0)| over N<=10 coefficients = {mp.nstr(diff, 3)}')

# ---------------------------------------------------------------- (d) pointwise, direct Legendre functions of complex degree
mp.mp.dps = 30
def Ghat_point(B, l0, r, eta):
    s = 0
    cache = {}
    for N, d in enumerate(B):
        for m, c in d.items():
            if m not in cache:
                cache[m] = mp.legenp(l0 + m, 0, eta, type=3)
            s += c * r**N * cache[m]
    return s
pts = [(mp.mpf('0.2'), mp.cosh(mp.mpf('0.6'))), (mp.mpf('0.3'), mp.cosh(mp.mpf('0.4'))), (mp.mpf('0.1'), mp.cosh(mp.mpf('1.5')))]
for nu in (mp.mpf('0.7'), mp.mpc('2.3', '0.2')):
    lv = mp.mpf(-0.5) + 1j * nu
    worst = 0
    for bb in range(-4, 5):
        B1 = recursion(lv, bb, 1, 45, conv_mp); Bm = recursion(lv - 1, bb, 0, 45, conv_mp); Bp = recursion(lv + 1, bb, 0, 45, conv_mp)
        for (r, eta) in pts:
            L = 2 * eta * (1 + r**2) * Ghat_point(B1, lv, r, eta)
            R = (1 + 1j / (2 * nu)) * Ghat_point(Bm, lv - 1, r, eta) + (1 - 1j / (2 * nu)) * Ghat_point(Bp, lv + 1, r, eta)
            worst = max(worst, abs(L - R) / abs(L))
    log(f'(d) pointwise (r, eta) in 3 points, b = -4..4, N <= 45, nu = {mp.nstr(nu, 4)}: max |LHS-RHS|/|LHS| = {mp.nstr(worst, 3)}')
log(f'done [{time.time()-t0:.0f}s]')
