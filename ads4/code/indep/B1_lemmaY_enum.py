"""Independent enumeration for Lemma Y (R4ARC_HEFT §2.3).  Written from the statement only.

Structures: monomials in the embedding invariants of four spin-2 legs,
    Z_i.Z_j (i<j), Z_i.P_j (i!=j), P_i.P_j (i<j)
(Z_i.Z_i = Z_i.P_i = P_i.P_i = 0), each Z_i of total degree exactly 2, times
at most one epsilon(V1..V5) of five DISTINCT vectors from {P_i, Z_i}
(a Z_i inside the epsilon uses one of its two legs).

For each structure:
    deg_i   = P_i-degree (monomial + epsilon slots)
    Delta'_i = 3 + deg_i ,  sigma = (Delta'_1 + Delta'_3 - Delta'_2 - Delta'_4)/2
    m_pol   = #{Z1.Z3, Z1.P3, Z3.P1, P1.P3}  (+1/2 if the epsilon holds exactly 3 of {P1,Z1,P3,Z3})
    eps with all 4 of {P1,Z1,P3,Z3} vanishes identically (P1-P3 = -t1 Z1 - t3 Z3) -> counted separately.
Lemma Y claims max(0, sigma) <= m_pol for all.

Several counting conventions for the P.P factors are tabulated, to compare with the quoted 997 836.
"""
import itertools, sys, os
from fractions import Fraction
from collections import Counter

OUT = os.path.join(os.path.dirname(__file__), 'cache', 'B1_lemmaY_enum.txt')
lines = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); lines.append(s)

legs = [1, 2, 3, 4]
ZZ = [(i, j) for i in legs for j in legs if i < j]
ZP = [(i, j) for i in legs for j in legs if i != j]          # Z_i . P_j
PP = [(i, j) for i in legs for j in legs if i < j]

def zstructures(zdeg):
    """all multisets of ZZ and ZP factors with Z_i degree = zdeg[i]."""
    out = []
    # choose ZZ multiplicities first (each <= 2)
    def rec_zz(k, cur, rem):
        if k == len(ZZ):
            rec_zp(0, cur, [], rem)
            return
        i, j = ZZ[k]
        for m in range(0, min(rem[i], rem[j]) + 1):
            r2 = dict(rem); r2[i] -= m; r2[j] -= m
            rec_zz(k + 1, cur + [((i, j), m)] if m else cur, r2)
    def rec_zp(k, zz, zp, rem):
        if k == len(ZP):
            if all(rem[i] == 0 for i in legs):
                out.append((tuple(zz), tuple(zp)))
            return
        i, j = ZP[k]
        # remaining ZP slots of leg i must be filled by later ZP entries with same i
        for m in range(0, rem[i] + 1):
            r2 = dict(rem); r2[i] -= m
            rec_zp(k + 1, zz, zp + [((i, j), m)] if m else zp, r2)
    rec_zz(0, [], dict(zdeg))
    return out

VEC = [('P', i) for i in legs] + [('Z', i) for i in legs]

def analyse(zz, zp, eps, pp):
    deg = {i: 0 for i in legs}
    for (i, j), m in zp: deg[j] += m
    for (i, j), m in pp: deg[i] += m; deg[j] += m
    if eps:
        for t, i in eps:
            if t == 'P': deg[i] += 1
    sigma = Fraction(deg[1] + deg[3] - deg[2] - deg[4], 2)
    m = 0
    for (i, j), mm in zz:
        if (i, j) == (1, 3): m += mm
    for (i, j), mm in zp:
        if (i, j) in ((1, 3), (3, 1)): m += mm
    for (i, j), mm in pp:
        if (i, j) == (1, 3): m += mm
    m = Fraction(m)
    k13 = 0
    if eps:
        k13 = sum(1 for t, i in eps if i in (1, 3))
        if k13 == 3: m += Fraction(1, 2)
    return sigma, m, k13

def pp_choices(rule):
    if rule == 'each<=2':
        for ms in itertools.product(range(3), repeat=6):
            yield tuple((PP[k], ms[k]) for k in range(6) if ms[k])
    elif rule == 'total<=2':
        for ms in itertools.product(range(3), repeat=6):
            if sum(ms) <= 2:
                yield tuple((PP[k], ms[k]) for k in range(6) if ms[k])
    elif rule == 'none':
        yield ()

# Z-structures
even_Z = zstructures({i: 2 for i in legs})
odd_Z = []   # (eps, zz, zp)
for eps in itertools.combinations(VEC, 5):
    zdeg = {i: 2 - sum(1 for t, j in eps if t == 'Z' and j == i) for i in legs}
    for zz, zp in zstructures(zdeg):
        odd_Z.append((eps, zz, zp))
say('parity-even Z-leg structures:', len(even_Z))
say('one-epsilon Z-leg structures (5 distinct of 8 vectors, 56 placements):', len(odd_Z))
k13c = Counter(sum(1 for t, i in e if i in (1, 3)) for e, _, _ in odd_Z)
say('   by k13 = #slots from {P1,Z1,P3,Z3}:', dict(sorted(k13c.items())))
say('   of which identically zero (k13 = 4):', k13c[4])

for rule in ['none', 'total<=2']:
    npp = sum(1 for _ in pp_choices(rule))
    worst = None; nbad = 0; tot = 0; tot_nz = 0
    wc = Counter()
    for zz, zp in even_Z:
        for pp in pp_choices(rule):
            s, m, _ = analyse(zz, zp, None, pp)
            mar = max(0, s) - m
            tot += 1; tot_nz += 1
            worst = mar if worst is None else max(worst, mar)
            nbad += mar > 0; wc[('even', mar)] += 1
    worst_e = worst
    worst = None
    for eps, zz, zp in odd_Z:
        k13 = sum(1 for t, i in eps if i in (1, 3))
        for pp in pp_choices(rule):
            s, m, _ = analyse(zz, zp, eps, pp)
            tot += 1
            if k13 == 4:
                continue
            tot_nz += 1
            mar = max(0, s) - m
            worst = mar if worst is None else max(worst, mar)
            nbad += mar > 0
    say(f'P.P rule {rule!r:12}: #P.P choices {npp:4d}; total structures {tot:8d} '
        f'(non-vanishing {tot_nz:8d}); worst margin max(0,sigma)-m: even {worst_e}, eps {worst}; violations {nbad}')

# 'each<=2' (3^6 = 729 P.P monomials): sigma and m depend on the P.P content only through
# f = mult(P1.P3) and f' = mult(P2.P4); the other 4 pairs contribute nothing (3^4 = 81 copies).
worst = None; nbad = 0; tot = 0; tot_nz = 0
for eps, zz, zp in [(None, a, b) for a, b in even_Z] + odd_Z:
    k13 = sum(1 for t, i in eps if i in (1, 3)) if eps else 0
    for f in range(3):
        for fp in range(3):
            pp = tuple(x for x in (((1, 3), f), ((2, 4), fp)) if x[1])
            s, m, _ = analyse(zz, zp, eps, pp)
            tot += 81
            if eps and k13 == 4: continue
            tot_nz += 81
            mar = max(0, s) - m
            worst = mar if worst is None else max(worst, mar); nbad += 81 * (mar > 0)
say(f"P.P rule 'each<=2'  : #P.P choices  729; total structures {tot:8d} (non-vanishing {tot_nz:8d}); worst margin {worst}; violations {nbad}")
# finer: margin distribution for even, rule 'none'
dist = Counter()
for zz, zp in even_Z:
    s, m, _ = analyse(zz, zp, None, ())
    dist[max(0, s) - m] += 1
say('even, no P.P: distribution of max(0,sigma)-m:', dict(sorted(dist.items())))
# check the exact identity sigma = m - (z + e' + f') + (k13-k24)/2  on all structures (no P.P and each<=2 P.P)
bad = 0
for eps, zz, zp in [(None, a, b) for a, b in even_Z] + odd_Z:
    for pp in pp_choices('total<=2'):
        s, m, k13 = analyse(zz, zp, eps, pp)
        z = sum(mm for (i, j), mm in zz if (i, j) == (2, 4))
        ep = sum(mm for (i, j), mm in zp if (i, j) in ((2, 4), (4, 2)))
        fp = sum(mm for (i, j), mm in pp if (i, j) == (2, 4))
        g = sum(mm for (i, j), mm in zz if (i, j) == (1, 3))
        e = sum(mm for (i, j), mm in zp if (i, j) in ((1, 3), (3, 1)))
        f = sum(mm for (i, j), mm in pp if (i, j) == (1, 3))
        kk = k13 if eps else 0
        k24 = (5 - kk) if eps else 0
        pred = Fraction(g + e + f - z - ep - fp) + Fraction(kk - k24, 2)
        bad += (pred != s)
say('closed form sigma = (g+e+f) - (z+e\'+f\') + (k13-k24)/2 : mismatches', bad)
open(OUT, 'w').write('\n'.join(lines) + '\n')
