"""Independent check of SM Theorem S18.9 (finite graph theorem) and the telescoping inverse (T9).

Random nonnegative weights w_ij on pairs of atoms define D via (T9); we check (a) that (T8) recovers w from D,
(b) that the min-cut of the chord-arrangement graph equals the non-crossing-matching minimum F_D(A) for every
union A of atoms. Control: the boundary-only Bell graph (capacities w_ij directly between atoms) must fail for
some disconnected unions.
"""
import itertools
import math
import random
from collections import deque

random.seed(7)


def nc_matchings(pts):
    if not pts:
        yield []
        return
    a = pts[0]
    for k in range(1, len(pts), 2):
        b = pts[k]
        inside, outside = pts[1:k], pts[k + 1:]
        for m1 in nc_matchings(inside):
            for m2 in nc_matchings(outside):
                yield [(a, b)] + m1 + m2


def maxflow(cap, s, t):
    # cap: dict of dict
    flow = 0.0
    res = {u: dict(v) for u, v in cap.items()}
    for u in list(res):
        for v in list(res[u]):
            res.setdefault(v, {}).setdefault(u, 0.0)
    while True:
        par = {s: None}
        q = deque([s])
        while q and t not in par:
            u = q.popleft()
            for v, c in res[u].items():
                if c > 1e-15 and v not in par:
                    par[v] = u
                    q.append(v)
        if t not in par:
            return flow
        # bottleneck
        v, bott = t, float('inf')
        while par[v] is not None:
            bott = min(bott, res[par[v]][v])
            v = par[v]
        v = t
        while par[v] is not None:
            u = par[v]
            res[u][v] -= bott
            res[v][u] += bott
            v = u
        flow += bott


def trial(n):
    th = [2 * math.pi * a / n + random.uniform(-0.2, 0.2) * 2 * math.pi / n for a in range(n)]
    P = [(math.cos(x), math.sin(x)) for x in th]
    # random weights (some zero), atoms i = arc (p_i, p_{i+1})
    w = {}
    for i in range(n):
        for j in range(i + 1, n):
            w[(i, j)] = 0.0 if random.random() < 0.25 else random.random()
    W = lambda i, j: w[(min(i, j), max(i, j))]

    def I(a, b):  # atoms on oriented arc p_a -> p_b
        out, k = [], a
        while k != b:
            out.append(k)
            k = (k + 1) % n
        return out
    D = [[0.0] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            if a != b:
                Ia = I(a, b)
                D[a][b] = sum(W(i, j) for i in Ia for j in range(n) if j not in Ia)
    # (a) T8 recovers w
    err8 = 0.0
    for i in range(n):
        for j in range(n):
            if i != j:
                ww = 0.5 * (D[i][j] + D[(i + 1) % n][(j + 1) % n] - D[i][(j + 1) % n] - D[(i + 1) % n][j])
                err8 = max(err8, abs(ww - W(i, j)))
    chords = [(a, b) for a in range(n) for b in range(a + 1, n)]

    def sgn(q):
        out = []
        for (a, b) in chords:
            ax, ay = P[a]
            bx, by = P[b]
            out.append(1 if (bx - ax) * (q[1] - ay) - (by - ay) * (q[0] - ax) > 0 else 0)
        return tuple(out)
    cap = {}
    term = {}
    for i in range(n):
        for j in range(i + 1, n):
            if W(i, j) == 0:
                continue
            # test chord from a point of atom i to a point of atom j
            def apt(k):
                t0, t1 = th[k], th[(k + 1) % n]
                if t1 < t0:
                    t1 += 2 * math.pi
                x = random.uniform(t0 + 0.1 * (t1 - t0), t1 - 0.1 * (t1 - t0))
                return (math.cos(x), math.sin(x))
            A, B = apt(i), apt(j)
            # crossing parameters with all chords
            ts = []
            for (a, b) in chords:
                ax, ay = P[a]
                bx, by = P[b]
                dx, dy = B[0] - A[0], B[1] - A[1]
                ex, ey = bx - ax, by - ay
                den = dx * ey - dy * ex
                if abs(den) < 1e-14:
                    continue
                s = ((ax - A[0]) * ey - (ay - A[1]) * ex) / den
                u = ((ax - A[0]) * dy - (ay - A[1]) * dx) / den
                if 1e-9 < s < 1 - 1e-9 and 1e-9 < u < 1 - 1e-9:
                    ts.append(s)
            ts = sorted(ts)
            grid = [0.0] + ts + [1.0]
            cells = []
            for k in range(len(grid) - 1):
                s = 0.5 * (grid[k] + grid[k + 1])
                cells.append(sgn((A[0] + s * (B[0] - A[0]), A[1] + s * (B[1] - A[1]))))
            term.setdefault(i, cells[0])
            term.setdefault(j, cells[-1])
            assert term[i] == cells[0] and term[j] == cells[-1]
            for k in range(len(cells) - 1):
                u, v = cells[k], cells[k + 1]
                cap.setdefault(u, {}).setdefault(v, 0.0)
                cap.setdefault(v, {}).setdefault(u, 0.0)
                cap[u][v] += W(i, j)
                cap[v][u] += W(i, j)
    maxerr, bellfail = 0.0, 0
    for r in range(1, n):
        for Aset in itertools.combinations(range(n), r):
            Aset = set(Aset)
            dA = [p for p in range(n) if ((p - 1) % n in Aset) != (p in Aset)]
            F = min(sum(D[a][b] for a, b in mt) for mt in nc_matchings(dA))
            g = {u: dict(v) for u, v in cap.items()}
            g['S'], g['T'] = {}, {}
            for i in range(n):
                if i in term:
                    if i in Aset:
                        g['S'][term[i]] = 1e9
                    else:
                        g.setdefault(term[i], {})['T'] = 1e9
            g.setdefault('T', {})
            mc = maxflow(g, 'S', 'T')
            maxerr = max(maxerr, abs(mc - F))
            bell = sum(W(i, j) for i in Aset for j in range(n) if j not in Aset)
            if bell < F - 1e-9:
                bellfail += 1
    return err8, maxerr, bellfail


if __name__ == '__main__':
    tot = 0
    for n in range(3, 9):
        for rep in range(4):
            e8, me, bf = trial(n)
            tot += 1
            print(f'n={n} rep={rep}: T8 inverse err={e8:.1e}  max|mincut-F_D|={me:.1e}  bell-graph failures={bf}')
    print('trials', tot)
