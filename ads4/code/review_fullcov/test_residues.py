"""PROOF.md §2 residue facts, checked against exact rationals / mpmath.
 (a) exact residue rho_{j,k} = lim_{h->1/p_j}(1-p_j h) U_k equals the closed formula
 (b) the partial-fraction decomposition U_k = U_k(inf) + sum_j rho_{j,k}/(1-p_j h) holds
 (c) |rho_{j,k}| <= (2/j) c^{-(s-1)} exp(-L(j,k)), L = 4-term lower bound (full_preflight.lower_exp)
 (d) |rho_{j,k}| increases in k to the limit and the limit is <= (2/j)c^{-(s-1)}
"""
import sys
from fractions import Fraction as Fr
import mpmath as mp
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'sliver_fullcov' / 'sliver_close'))
import os; os.environ.setdefault('SLIVER_COEFF', __import__('rc').COEFF)
from full_preflight import lower_exp

mp.mp.dps = 60

def rho_exact(j, k, c, s):
    p = lambda i: i * (i + s)
    pj = p(j); r = c - 1
    for i in range(k):
        if i != j:
            r *= (c * pj - p(i)) / Fr(pj - p(i))
    return r

def rho_formula(j, k, c, s):
    pj = j * (j + s)
    v = (-s + mp.sqrt(s * s + 4 * mp.mpf(c.numerator) / c.denominator * pj)) / 2
    return (2 * j + s) / mp.mpf(pj) * abs(mp.rf(-v, k) * mp.rf(v + s, k)) / (mp.factorial(j) * mp.factorial(k - 1 - j) * mp.rf(j + s, k))

def rho_limit(j, c, s):
    pj = j * (j + s)
    v = (-s + mp.sqrt(s * s + 4 * mp.mpf(c.numerator) / c.denominator * pj)) / 2
    return (2 * j + s) / mp.mpf(pj) * mp.gamma(j + s) * abs(mp.rgamma(-v)) / (mp.factorial(j) * mp.gamma(v + s))

worst_a = 0; worst_c = 0; worst_d = 0; nb = 0
for s in (1, 9):
    for c in (Fr(258, 1000), Fr(1, 2), Fr(3, 4), Fr(9, 10), Fr(899, 1000)):
        for j in (1, 2, 3, 5, 8, 13, 30, 77, 150):
            lim = rho_limit(j, c, s)
            bound_lim = mp.mpf(2) / j * (mp.mpf(c.numerator) / c.denominator) ** (-(s - 1))
            worst_d = max(worst_d, lim / bound_lim)
            prev = 0
            for k in sorted(set([j + 1, j + 2, j + 3, 2 * j + 1, 3 * j, 5 * j, 10 * j + 3, 40 * j, 400])):
                if k <= j: continue
                ex = rho_exact(j, k, c, s); exm = abs(mp.mpf(ex.numerator) / ex.denominator)
                fo = rho_formula(j, k, c, s)
                worst_a = max(worst_a, abs(exm - fo) / exm) if exm else max(worst_a, fo)
                assert exm >= prev * (1 - mp.mpf(10) ** -40), 'not increasing in k'
                prev = exm
                L = lower_exp(j, k, c, s)
                bd = bound_lim * mp.exp(-mp.mpf(L.numerator) / L.denominator)
                worst_c = max(worst_c, exm / bd); nb += 1
                assert exm <= bd * (1 + mp.mpf(10) ** -40), ('residue bound violated', s, c, j, k, exm, bd)
                assert exm <= lim * (1 + mp.mpf(10) ** -40)
print('max rel diff exact vs closed formula:', mp.nstr(worst_a, 5))
print('max |rho_jk| / [(2/j)c^-(s-1) e^-L(j,k)] over', nb, 'cases:', mp.nstr(worst_c, 8))
print('max |rho_inf| / [(2/j)c^-(s-1)]:', mp.nstr(worst_d, 8))

# (b) partial fractions at small k, exact rationals
for s in (1, 9):
    for c in (Fr(3, 10), Fr(7, 8)):
        for k in (1, 2, 5, 9):
            for h in (Fr(1, 1000), Fr(3, 77), Fr(-2, 5)):
                U = Fr(1)
                for jj in range(k): U *= (c - jj * (jj + s) * h) / (1 - jj * (jj + s) * h)
                Uinf = c if k >= 1 else Fr(1)
                PF = Uinf + sum(rho_exact(jj, k, c, s) / (1 - jj * (jj + s) * h) for jj in range(1, k))
                assert PF == U, (s, c, k, h)
print('partial-fraction decomposition exact: PASS')
