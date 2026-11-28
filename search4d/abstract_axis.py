"""Abstract test: entropy vectors from the 'matching formula with a free axis'.

Setting (SO(2)-reduced AdS4 problem): marked boundary edges e_1<...<e_m on the
ray rho>0, the axis is a free arc next to the corner rho=0.  Data:
  L[i,j]  = area of the extremal tube joining edges i,j
  c[i]    = area of the extremal cap from edge i to the axis
S(X) = min over non-crossing configurations of tubes and caps, where caps may
not lie under a tube, and any number of caps may end on the axis.

Questions tested here:
  Q1 (disc sanity, 2+1 analogue): caps forbidden, L Kalmanson -> in RT cone?
  Q2 L,c with augmented circle {0,e_1..e_m} circular-decomposable -> in cone?
  Q3 L,c obeying only the SSA/WM-derived linear conditions K1,K2,K4
     (K0: c_i+c_j>=L_ij allowed to FAIL) -> in cone?
  Q4 random L,c with NO conditions (control: violations should appear,
     showing the test has teeth).
Run: python3 -I abstract_axis.py [nsamples]
"""
import json, sys, itertools, os
import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
FJ = json.load(open(os.path.join(HERE, '..', 'strips', 'data', 'facets_n5_all372.json')))
SUBS = FJ['subset_order']
F = np.array(FJ['facets'], dtype=float)
ORB = FJ['facet_orbits']
ORBIT_OF = np.zeros(len(F), int)
for k, o in enumerate(ORB):
    ORBIT_OF[o] = k
SUBSETS = [frozenset('ABCDE'.index(ch) for ch in s) for s in SUBS]


def balanced_check():
    labels = range(6)  # 5 = purifier
    bad = []
    for fi, q in enumerate(F):
        for x, y in itertools.combinations(labels, 2):
            tot = sum(q[i] * ((x in S) != (y in S)) for i, S in enumerate(SUBSETS))
            if abs(tot) > 1e-12:
                bad.append(ORBIT_OF[fi]); break
    return sorted(set(bad))


INF = 1e300


def entropy(active, L, c):
    """active: sorted list of edge indices; returns min configuration area."""
    k = len(active)
    if k == 0:
        return 0.0
    P = active
    # M[a][b]: min noncrossing perfect matching of P[a..b-1] (half-open)
    M = [[INF] * (k + 1) for _ in range(k + 1)]
    for a in range(k + 1):
        M[a][a] = 0.0
    for length in range(2, k + 1, 2):
        for a in range(0, k - length + 1):
            b = a + length
            best = INF
            for cc in range(a + 1, b, 2):
                val = L[P[a], P[cc]] + M[a + 1][cc] + M[cc + 1][b]
                if val < best:
                    best = val
            M[a][b] = best
    # caps: g[j] = min cost of P[0..j-1] with a cap at P[j-1] (g[0]=0 base)
    g = [INF] * (k + 1)
    g[0] = 0.0
    for j in range(1, k + 1):
        if c[P[j - 1]] >= INF:
            continue
        best = INF
        for i in range(0, j):
            if g[i] < INF and M[i][j - 1] < INF:
                v = g[i] + M[i][j - 1]
                if v < best:
                    best = v
        g[j] = best + c[P[j - 1]] if best < INF else INF
    best = INF
    for i in range(0, k + 1):
        if g[i] < INF and M[i][k] < INF:
            best = min(best, g[i] + M[i][k])
    return best


def entropy_vector(labels, L, c):
    """labels: party label (0..4, 5=O) of intervals I_0..I_m (I_m=[e_m,inf)=O)."""
    m = len(labels) - 1
    vec = np.empty(31)
    for si, S in enumerate(SUBSETS):
        inside = [lab in S for lab in labels]
        active = [j for j in range(m) if inside[j] != inside[j + 1]]
        vec[si] = entropy(active, L, c)
    return vec


def random_labels(m, rng):
    while True:
        lab = [int(rng.integers(0, 6)) for _ in range(m)] + [5]
        if any(lab[i] == lab[i + 1] for i in range(m)):
            continue
        if set(range(5)) <= set(lab):
            return lab


def circ_decomp(m, rng, sparse=0.3, nocaps=False):
    """Points on circle: index 0 = axis, 1..m edges.  Random split weights."""
    n = m + 1
    D = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):  # arc {a..b-1}
            if rng.random() < sparse:
                w = rng.exponential()
                ins = np.zeros(n, bool); ins[a:b] = True
                D += w * (ins[:, None] != ins[None, :])
    L = D[1:, 1:].copy()
    c = D[0, 1:].copy()
    if nocaps:
        c[:] = INF
    return L, c


def kalmanson_lp_sample(m, rng, k0_free=True):
    """Random vertex of {L,c >=0, <=1 box, K1,K2,K4} via LP with random objective."""
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    pid = {p: t for t, p in enumerate(pairs)}
    nL = len(pairs); nv = nL + m
    def Lx(i, j):
        e = np.zeros(nv); e[pid[(min(i, j), max(i, j))]] = 1; return e
    def cx(i):
        e = np.zeros(nv); e[nL + i] = 1; return e
    A = []  # rows r: r.x >= 0
    for i, j, k in itertools.combinations(range(m), 3):
        A.append(Lx(i, k) + cx(j) - Lx(i, j) - cx(k))  # K1
        A.append(Lx(i, k) + cx(j) - cx(i) - Lx(j, k))  # K2
    for i, j, k, l in itertools.combinations(range(m), 4):
        A.append(Lx(i, k) + Lx(j, l) - Lx(i, j) - Lx(k, l))
        A.append(Lx(i, k) + Lx(j, l) - Lx(i, l) - Lx(j, k))
    if not k0_free:
        for i, j in pairs:
            A.append(cx(i) + cx(j) - Lx(i, j))
    A = np.array(A)
    obj = rng.normal(size=nv)
    res = linprog(obj, A_ub=-A, b_ub=np.zeros(len(A)), bounds=[(0, 1)] * nv, method='highs')
    x = res.x
    L = np.zeros((m, m))
    for (i, j), t in pid.items():
        L[i, j] = L[j, i] = x[t]
    c = x[nL:].copy()
    return L, c


def margins(vec):
    return F @ vec


def run(mode, nsamp, rng, mrange=(5, 10), nlab=30):
    worst = np.full(len(ORB), np.inf)
    worst_case = [None] * len(ORB)
    k0viol = 0
    for s in range(nsamp):
        m = int(rng.integers(mrange[0] + (1 if mode == 'disc' else 0), mrange[1] + 1))
        if mode == 'disc':
            L, c = circ_decomp(m, rng, nocaps=True)
        elif mode == 'circ':
            L, c = circ_decomp(m, rng)
        elif mode == 'K124':
            L, c = kalmanson_lp_sample(m, rng, k0_free=True)
            if any(L[i, j] > c[i] + c[j] + 1e-9 for i in range(m) for j in range(i + 1, m)):
                k0viol += 1
        elif mode == 'inflate':
            # physical constraint is only on L' = min(L, c_i+c_j): take L' circular-decomposable
            # (sparse split system so that L' = c_i+c_j happens) and inflate raw L there.
            L, c = circ_decomp(m, rng, sparse=0.12)
            for i in range(m):
                for j in range(i + 1, m):
                    if L[i, j] >= c[i] + c[j] - 1e-12:
                        L[i, j] = L[j, i] = c[i] + c[j] + (rng.exponential() if rng.random() < 0.7 else 1e6)
                        k0viol += 1
        elif mode == 'mix':
            L1, c1 = kalmanson_lp_sample(m, rng, k0_free=True)
            L2, c2 = circ_decomp(m, rng)
            w = rng.random()
            L = w * L1 + (1 - w) * L2 / max(1e-9, L2.max()); c = w * c1 + (1 - w) * c2 / max(1e-9, L2.max())
        elif mode == 'random':
            L = rng.random((m, m)); L = (L + L.T) / 2; c = rng.random(m)
        for _ in range(nlab):
            lab = random_labels(m, rng)
            while mode == 'disc' and lab[0] != 5:
                # no axis: the corner interval must also be purifier (line, not half-line)
                lab = random_labels(m, rng)
            vec = entropy_vector(lab, L, c)
            if not np.all(np.isfinite(vec)) or np.any(vec >= INF / 2):
                continue
            mg = margins(vec)
            for o in range(len(ORB)):
                if o == 2:
                    continue  # SA: not balanced, needs endpoint terms
                v = mg[ORB[o]].min()
                if v < worst[o]:
                    worst[o] = v; worst_case[o] = (m, lab, ORB[o][int(np.argmin(mg[ORB[o]]))])
    return worst, worst_case, k0viol


if __name__ == '__main__':
    nsamp = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    rng = np.random.default_rng(1)
    print('orbits needing endpoint terms (unbalanced):', balanced_check())
    modes = sys.argv[2].split(',') if len(sys.argv) > 2 else ['disc', 'circ', 'K124', 'inflate', 'mix', 'random']
    for mode in modes:
        w, wc, k0 = run(mode, nsamp, rng)
        print(f'{mode:7s} min margin per orbit (SA skipped):',
              ' '.join(f'{x:+.3e}' if np.isfinite(x) else '   n/a   ' for x in w),
              f' [K0-violating samples: {k0}]' if mode in ('K124', 'inflate') else '', flush=True)
