"""Count independent LOCAL five-point superamplitude contacts in 4D N=8 supergravity by mass dimension (EFK-style test).
MHV sector: M5 = delta^16(Q) f, with f totally S5-symmetric and of little-group weight t_i^{-4} for every leg, and the component
P := M(1-2-3+4+5+) = <12>^8 f must be a polynomial in spinor brackets of mass dimension D (= number of derivatives of the operator:
R^5 <-> D = 10, D^2 R^5 <-> 12, ...).  For n = 5 the anti-MHV sector is the parity conjugate (same count).
Method: enumerate all bracket monomials with weights (w1..w5) = (+4,+4,-4,-4,-4) and dimension D, evaluate on random complex
momentum-conserving kinematics, find the function space (rank), and impose f(sigma kin) = f(kin) for sigma in {(12), (12345)}.
Output: number of independent local supersymmetric five-point contacts at each D."""
import itertools, sys
import numpy as np
rng = np.random.default_rng(7)
PAIRS = list(itertools.combinations(range(5), 2))
W = np.array([4, 4, -4, -4, -4])

def kin():
    lam = rng.normal(size=(5, 2)) + 1j*rng.normal(size=(5, 2))
    lt = rng.normal(size=(5, 2)) + 1j*rng.normal(size=(5, 2))
    # solve lt3, lt4 (indices 3,4) from sum_i lam_i (x) lt_i = 0
    S = sum(np.outer(lam[i], lt[i]) for i in range(3))
    A = np.column_stack([lam[3], lam[4]])          # 2x2
    X = np.linalg.solve(A, -S)                     # rows: lt3, lt4
    lt[3], lt[4] = X[0], X[1]
    assert np.abs(sum(np.outer(lam[i], lt[i]) for i in range(5))).max() < 1e-9
    return lam, lt

def brackets(lam, lt):
    ang = np.array([lam[i][0]*lam[j][1] - lam[i][1]*lam[j][0] for (i, j) in PAIRS])
    sq = np.array([lt[i][0]*lt[j][1] - lt[i][1]*lt[j][0] for (i, j) in PAIRS])
    return ang, sq

def multigraphs(degrees, nedges):
    """all multisets of pairs with given vertex degrees (exponent vectors over PAIRS)."""
    out = []
    deg = list(degrees)
    def rec(k, rem, cur):
        if k == len(PAIRS):
            if all(x == 0 for x in rem): out.append(tuple(cur))
            return
        i, j = PAIRS[k]
        m = min(rem[i], rem[j])
        # remaining pairs must cover: prune if vertex i has no later pairs
        for e in range(m, -1, -1):
            rem[i] -= e; rem[j] -= e; cur.append(e)
            ok = True
            # vertex v with remaining degree but no later pair containing v -> prune
            for v in range(5):
                if rem[v] > 0 and not any(v in PAIRS[q] for q in range(k + 1, len(PAIRS))):
                    ok = False; break
            if ok: rec(k + 1, rem, cur)
            cur.pop(); rem[i] += e; rem[j] += e
    if any(d < 0 for d in deg) or sum(deg) != 2*nedges: return []
    rec(0, deg[:], [])
    return out

def monomials(D):
    Na = (D - 2)//2; Nb = Na + 2
    mons = []
    for a in itertools.combinations_with_replacement(range(10), Na):
        av = np.bincount(a, minlength=10)
        Adeg = np.zeros(5, int)
        for k, (i, j) in enumerate(PAIRS):
            Adeg[i] += av[k]; Adeg[j] += av[k]
        Bdeg = Adeg - W
        for bv in multigraphs(list(Bdeg), Nb):
            mons.append((tuple(av), bv))
    return mons

def evalmons(mons, ang, sq):
    A = np.array([m[0] for m in mons]); B = np.array([m[1] for m in mons])
    return np.prod(ang[None, :]**A, axis=1)*np.prod(sq[None, :]**B, axis=1)

def count(D, npts=None):
    mons = monomials(D)
    nm = len(mons)
    npts = npts or max(3*nm, 60)
    perms = [(1, 0, 2, 3, 4), (1, 2, 3, 4, 0)]
    rows, crows = [], []
    for _ in range(npts):
        lam, lt = kin()
        ang, sq = brackets(lam, lt)
        e0 = evalmons(mons, ang, sq)/ang[0]**8                    # f values
        rows.append(e0)
        for p in perms:
            ang2, sq2 = brackets(lam[list(p)], lt[list(p)])
            crows.append(evalmons(mons, ang2, sq2)/ang2[0]**8 - e0)
    E = np.array(rows); C = np.array(crows)
    sc = np.abs(E).max(axis=0); sc[sc == 0] = 1
    E /= sc; C /= sc
    u, s, vh = np.linalg.svd(E, full_matrices=False)
    r = int(np.sum(s > s[0]*1e-9))
    Bf = vh[:r].conj().T                     # coefficient vectors spanning the function space (mod relations)
    Cr = C @ Bf
    s2 = np.linalg.svd(Cr, compute_uv=False)
    r2 = int(np.sum(s2 > s2[0]*1e-8)) if len(s2) and s2[0] > 0 else 0
    return nm, r, r - r2, s2[-3:] if len(s2) >= 3 else s2

PR = 2147483647
def kin_mod():
    while True:
        lam = [[int(rng.integers(1, PR)) for _ in range(2)] for _ in range(5)]
        lt = [[int(rng.integers(1, PR)) for _ in range(2)] for _ in range(5)]
        S = [[sum(lam[i][a]*lt[i][b] for i in range(3)) % PR for b in range(2)] for a in range(2)]
        det = (lam[3][0]*lam[4][1] - lam[3][1]*lam[4][0]) % PR
        if det == 0: continue
        inv = pow(det, PR - 2, PR)
        # solve [lam3 lam4] X = -S  (columns lam3, lam4): X = A^-1 (-S)
        Ainv = [[lam[4][1]*inv % PR, -lam[4][0]*inv % PR], [-lam[3][1]*inv % PR, lam[3][0]*inv % PR]]
        for b in range(2):
            x3 = (-(Ainv[0][0]*S[0][b] + Ainv[0][1]*S[1][b])) % PR
            x4 = (-(Ainv[1][0]*S[0][b] + Ainv[1][1]*S[1][b])) % PR
            lt[3][b], lt[4][b] = x3, x4
        for a in range(2):
            for b in range(2):
                assert sum(lam[i][a]*lt[i][b] for i in range(5)) % PR == 0
        return lam, lt

def br_mod(lam, lt):
    ang = [(lam[i][0]*lam[j][1] - lam[i][1]*lam[j][0]) % PR for (i, j) in PAIRS]
    sq = [(lt[i][0]*lt[j][1] - lt[i][1]*lt[j][0]) % PR for (i, j) in PAIRS]
    return ang, sq

def ev_mod(mons, ang, sq):
    inv8 = pow(pow(ang[0], 8, PR), PR - 2, PR)
    out = []
    for (av, bv) in mons:
        v = 1
        for k in range(10):
            if av[k]: v = v*pow(ang[k], int(av[k]), PR) % PR
            if bv[k]: v = v*pow(sq[k], int(bv[k]), PR) % PR
        out.append(v*inv8 % PR)
    return out

def count_mod(D):
    import flint
    mons = monomials(D); nm = len(mons)
    npts = nm + 40
    perms = [(1, 0, 2, 3, 4), (1, 2, 3, 4, 0)]
    E, C = [], []
    for _ in range(npts):
        lam, lt = kin_mod()
        e0 = ev_mod(mons, *br_mod(lam, lt)); E.append(e0)
        for p in perms:
            e1 = ev_mod(mons, *br_mod([lam[i] for i in p], [lt[i] for i in p]))
            C.append([(x - y) % PR for x, y in zip(e1, e0)])
    rE = flint.nmod_mat(len(E), nm, [x for row in E for x in row], PR).rank()
    rC = flint.nmod_mat(len(C), nm, [x for row in C for x in row], PR).rank()
    # nsym = nullity(C) - nullity(E) = (nm - rC) - (nm - rE) = rE - rC
    return nm, rE, rE - rC

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'mod':
    for D in [int(x) for x in sys.argv[2:]]:
        nm, r, nsym = count_mod(D)
        print('[exact mod p] D=%2d: monomials %6d, independent (--+++) polynomials %5d, S5-symmetric local N=8 superamplitudes (MHV sector): %d' % (D, nm, r, nsym))
        sys.stdout.flush()
    sys.exit(0)

if __name__ == '__main__':
    for D in [int(x) for x in sys.argv[1:]] or [10, 12, 14, 16, 18]:
        nm, r, nsym, tail = count(D)
        print('D=%2d: monomials %6d, independent polynomials (helicity --+++) %5d, S5-symmetric local superamplitudes: %d   (smallest sv %s)' % (
            D, nm, r, nsym, np.array2string(np.abs(tail), precision=2)))
        sys.stdout.flush()
