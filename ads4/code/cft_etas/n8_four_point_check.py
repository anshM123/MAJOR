"""Validation of the counting method on four points (known answer): local N=8 superamplitudes delta^16(Q) f with f S4-symmetric,
<12>^8 f polynomial of dimension D.  Expected: D=8 (R^4) 1; D=10 (D^2R^4) 0; D=12 (D^4R^4) 1; D=14 (D^6R^4) 1; D=16: 1; D=18: 1; D=20: 2."""
import itertools, sys
import numpy as np, flint
rng = np.random.default_rng(3)
PR = 2147483647
N = 4
PAIRS = list(itertools.combinations(range(N), 2))
W = [4, 4, -4, -4]
def kin_mod():
    while True:
        lam = [[int(rng.integers(1, PR)) for _ in range(2)] for _ in range(N)]
        lt = [[int(rng.integers(1, PR)) for _ in range(2)] for _ in range(N)]
        S = [[sum(lam[i][a]*lt[i][b] for i in range(N - 2)) % PR for b in range(2)] for a in range(2)]
        det = (lam[N-2][0]*lam[N-1][1] - lam[N-2][1]*lam[N-1][0]) % PR
        if det == 0: continue
        inv = pow(det, PR - 2, PR)
        Ainv = [[lam[N-1][1]*inv % PR, -lam[N-1][0]*inv % PR], [-lam[N-2][1]*inv % PR, lam[N-2][0]*inv % PR]]
        for b in range(2):
            lt[N-2][b] = (-(Ainv[0][0]*S[0][b] + Ainv[0][1]*S[1][b])) % PR
            lt[N-1][b] = (-(Ainv[1][0]*S[0][b] + Ainv[1][1]*S[1][b])) % PR
        return lam, lt
def br(lam, lt):
    return ([(lam[i][0]*lam[j][1] - lam[i][1]*lam[j][0]) % PR for (i, j) in PAIRS],
            [(lt[i][0]*lt[j][1] - lt[i][1]*lt[j][0]) % PR for (i, j) in PAIRS])
def monomials(D):
    out = []
    np_ = len(PAIRS)
    for tot_a in range(D + 1):
        for a in itertools.combinations_with_replacement(range(np_), tot_a):
            av = [a.count(k) for k in range(np_)]
            for b in itertools.combinations_with_replacement(range(np_), D - tot_a):
                bv = [b.count(k) for k in range(np_)]
                ok = all(sum(av[k] - bv[k] for k, (i, j) in enumerate(PAIRS) if v in (i, j)) == W[v] for v in range(N))
                if ok: out.append((av, bv))
    return out
def ev(mons, ang, sq):
    inv8 = pow(pow(ang[0], 8, PR), PR - 2, PR)
    res = []
    for av, bv in mons:
        v = 1
        for k in range(len(PAIRS)):
            if av[k]: v = v*pow(ang[k], av[k], PR) % PR
            if bv[k]: v = v*pow(sq[k], bv[k], PR) % PR
        res.append(v*inv8 % PR)
    return res
for D in (8, 10, 12, 14, 16, 18, 20):
    mons = monomials(D); nm = len(mons)
    E, C = [], []
    perms = [(1, 0, 2, 3), (1, 2, 3, 0)]
    for _ in range(nm + 30):
        lam, lt = kin_mod(); e0 = ev(mons, *br(lam, lt)); E.append(e0)
        for p in perms:
            e1 = ev(mons, *br([lam[i] for i in p], [lt[i] for i in p])); C.append([(x - y) % PR for x, y in zip(e1, e0)])
    rE = flint.nmod_mat(len(E), nm, [x for r in E for x in r], PR).rank()
    rC = flint.nmod_mat(len(C), nm, [x for r in C for x in r], PR).rank()
    print('4-pt D=%2d: monomials %4d, independent polys %3d, S4-symmetric local superamplitudes %d' % (D, nm, rE, rE - rC))
