"""Is 'subtract (1/2) I(alpha:beta) Bell pairs for cyclically adjacent pieces' cone-preserving for ARBITRARY
graph models? Random graphs: n boundary terminals (pieces on a circle, the last one may be the purifier),
k bulk vertices, random capacities. S_fine = min-cut. G' = S_fine - sum_{adjacent} (1/2) I b. Test SA and MMI
on G' for all party groupings of the pieces into <= N parties (adjacent pieces distinct)."""
import itertools, math, sys
import numpy as np
def mincut_all(W, n, k):
    V = n + k; best = {}
    for mask in range(1, 2**n - 1):
        T = [i for i in range(n) if mask >> i & 1]
        b = math.inf
        for bm in range(2**k):
            S = set(T) | {n + j for j in range(k) if bm >> j & 1}
            b = min(b, sum(W[i, j] for i in S for j in range(V) if j not in S))
        best[frozenset(T)] = b
    best[frozenset()] = 0.0; best[frozenset(range(n))] = 0.0
    return best
def run(trials, seed, n=5, k=4, N=3):
    rng = np.random.default_rng(seed); viol = 0; worst = 0
    for t in range(trials):
        V = n + k; W = np.zeros((V, V))
        for i in range(V):
            for j in range(i+1, V):
                if rng.random() < 0.6: W[i, j] = W[j, i] = rng.exponential(1.0)
        S = mincut_all(W, n, k)
        I = lambda a, b: S[frozenset({a})] + S[frozenset({b})] - S[frozenset({a, b})]
        tp = {p: 0.5*I((p-1) % n, p) for p in range(n)}      # point p sits between piece p-1 and piece p
        def G(X):   # X = set of pieces
            E = [p for p in range(n) if (((p-1) % n) in X) != (p in X)]
            return S[frozenset(X)] - sum(tp[p] for p in E)
        # fine-grained: each piece its own party (pure state on n pieces)
        pieces = range(n); sub = [frozenset(c) for r in range(1, n) for c in itertools.combinations(pieces, r)]
        bad = 0
        for X, Y in itertools.combinations(sub, 2):
            if X & Y or len(X | Y) == n: continue
            q = G(X) + G(Y) - G(X | Y); worst = min(worst, q); bad |= q < -1e-9
        for X, Y, Z in itertools.combinations(sub, 3):
            if X & Y or X & Z or Y & Z: continue
            q = G(X | Y) + G(Y | Z) + G(X | Z) - G(X) - G(Y) - G(Z) - G(X | Y | Z); worst = min(worst, q); bad |= q < -1e-9
        neg = min(G(X) for X in sub); bad |= neg < -1e-9; worst = min(worst, neg)
        viol += bad
    print(f"random graphs n={n} pieces, k={k} bulk vertices: {viol}/{trials} violate SA/MMI/positivity after "
          f"adjacent-MI subtraction; most negative value {worst:.3g}")
if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 300, 1)
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 300, 2, n=6, k=3)
