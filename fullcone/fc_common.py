"""Shared helpers: float min-cut (Edmonds-Karp), graph-model entropy vectors, facet checks, contraction maps."""
import itertools, collections, sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cones

def mincut(nodes, edges, src_set, snk_set):
    """edges: list of (u,v,w) undirected. min weight cut separating src_set from snk_set."""
    idx = {u: k for k, u in enumerate(nodes)}; n = len(nodes) + 2; S, T = n-2, n-1
    cap = collections.defaultdict(float); adj = collections.defaultdict(set)
    def add(a, b, w):
        cap[(a, b)] += w; cap[(b, a)] += w; adj[a].add(b); adj[b].add(a)
    for u, v, w in edges:
        if w > 0 and u != v: add(idx[u], idx[v], w)
    BIG = 1e18
    for u in src_set: add(S, idx[u], BIG)
    for u in snk_set: add(idx[u], T, BIG)
    flow = 0.0
    while True:
        par = {S: None}; q = collections.deque([S])
        while q and T not in par:
            a = q.popleft()
            for b in adj[a]:
                if b not in par and cap[(a, b)] > 1e-12:
                    par[b] = a; q.append(b)
        if T not in par: return flow
        f = BIG; b = T
        while par[b] is not None: f = min(f, cap[(par[b], b)]); b = par[b]
        b = T
        while par[b] is not None:
            a = par[b]; cap[(a, b)] -= f; cap[(b, a)] += f; b = a
        flow += f

def graph_entropy(nodes, edges, label, N):
    """label: node -> party (0..N-1, N=purifier) for terminals. returns dict mask->S"""
    out = {}
    for m in range(1, 2**N):
        src = [u for u in nodes if u in label and label[u] < N and (m >> label[u]) & 1]
        snk = [u for u in nodes if u in label and not (label[u] < N and (m >> label[u]) & 1)]
        out[m] = mincut(nodes, edges, src, snk)
    return out

_INST = {}
def facet_check(S, N, tol=1e-9):
    """returns (min over facets of value/scale, name of worst facet)"""
    if N not in _INST: _INST[N] = cones.instances(N)
    scale = max(1.0, max(abs(v) for v in S.values()))
    worst, wn = float("inf"), None
    for name, ineq in _INST[N]:
        v = cones.evaluate(ineq, S)/scale
        if v < worst: worst, wn = v, name
    return worst, wn

def contraction_exists(ineq, N):
    """ineq: tuple((mask,coeff)) over parties 0..N-1 (purifier N implicit). search a contraction map."""
    lhs = [(c, m) for m, c in ineq if c > 0]; rhs = [(-c, m) for m, c in ineq if c < 0]
    L, R = len(lhs), len(rhs); al = [c for c, _ in lhs]; be = [c for c, _ in rhs]
    occ = lambda terms, i: tuple(((m >> i) & 1) if i < N else 0 for _, m in terms)
    X = list(itertools.product((0, 1), repeat=L)); Y = list(itertools.product((0, 1), repeat=R))
    dX = lambda x, y: sum(a for a, s, t in zip(al, x, y) if s != t)
    dY = lambda x, y: sum(b for b, s, t in zip(be, x, y) if s != t)
    fixed = {}
    for i in range(N+1):
        x, y = occ(lhs, i), occ(rhs, i)
        if x in fixed and fixed[x] != y: return False
        fixed[x] = y
    for a, b in itertools.combinations(fixed, 2):
        if dY(fixed[a], fixed[b]) > dX(a, b): return False
    order = sorted(X, key=lambda x: (x not in fixed, min(dX(x, a) for a in fixed)))
    f = {}
    def rec(k):
        if k == len(order): return True
        x = order[k]
        for y in ([fixed[x]] if x in fixed else Y):
            if all(dY(y, f[z]) <= dX(x, z) for z in f):
                f[x] = y
                if rec(k+1): return True
                del f[x]
        return False
    return rec(0)
