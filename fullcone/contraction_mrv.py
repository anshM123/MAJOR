"""Independent contraction-map search (MRV + forward checking, numpy domains) for weighted inequalities.
Usage: python3 -I contraction_mrv.py  -> tests the N=5 orbit representatives of cones.py."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, itertools, cones

def search(ineq, N, limit=2_000_000):
    lhs = [(c, m) for m, c in ineq if c > 0]; rhs = [(-c, m) for m, c in ineq if c < 0]
    L, R = len(lhs), len(rhs)
    al = np.array([c for c, _ in lhs]); be = np.array([c for c, _ in rhs])
    X = np.array(list(itertools.product((0, 1), repeat=L))); Y = np.array(list(itertools.product((0, 1), repeat=R)))
    DX = (X[:, None, :] != X[None, :, :]) @ al
    DY = ((Y[:, None, :] != Y[None, :, :]) @ be).astype(np.int16)
    xi = {tuple(X[k]): k for k in range(len(X))}; yi = {tuple(Y[k]): k for k in range(len(Y))}
    dom = np.ones((len(X), len(Y)), bool)
    for i in range(N+1):
        x = tuple(((m >> i) & 1) if i < N else 0 for _, m in lhs); y = tuple(((m >> i) & 1) if i < N else 0 for _, m in rhs)
        d = np.zeros(len(Y), bool); d[yi[y]] = True; dom[xi[x]] &= d
    assigned = {}; nodes = [0]
    def propagate(k, v, dom):
        # after assigning x_k -> y_v, restrict others: DY[y, v] <= DX[x, k]
        ok = DY[:, v][None, :] <= DX[:, k][:, None]
        return dom & ok
    sys.setrecursionlimit(10000)
    def rec(dom):
        nodes[0] += 1
        if nodes[0] > limit: raise TimeoutError
        free = [k for k in range(len(X)) if k not in assigned]
        if not free: return True
        sizes = dom[free].sum(1)
        if sizes.min() == 0: return False
        k = free[int(np.argmin(sizes))]
        for v in np.flatnonzero(dom[k]):
            assigned[k] = v
            nd = propagate(k, v, dom); nd[k] = False; nd[k, v] = True
            if rec(nd): return True
            del assigned[k]
        return False
    try:
        return rec(dom), nodes[0]
    except TimeoutError:
        return None, nodes[0]

if __name__ == "__main__":
    reps = {}
    for name, ineq in cones.instances(5): reps.setdefault(name, ineq)
    for name in sys.argv[1:] or reps:
        r, n = search(reps[name], 5)
        print(name, "contraction map:", r, "nodes", n, flush=True)
