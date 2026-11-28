"""General 'disc orbit space with free arcs' entropy engine (copies model).

The SO(2)-orbit space of a Cauchy slice is a disc whose boundary circle is a
cyclic sequence of asymptotic arcs (carrying marked edges) and free arcs
(rotation axes: surfaces may end there at no cost).  An HRT configuration is a
non-crossing chord diagram: edge-edge chords (tubes), edge-axis chords (caps),
axis-axis chords between different axes (closed spheres, e.g. a horizon).
Model: each free arc is k copies; sub-arcs between copies take free membership
X; S(A) = min_X  min noncrossing perfect matching of the switch points.

points: list of ('e', i) or ('f', a) in cyclic order; arcs[k] lies between
points[k] and points[k+1]; label[k] in 0..5 for asymptotic arcs, None for free
sub-arcs.  Distances: D(e_i,e_j)=L[i,j], D(e_i,f_a)=C[a][i], D(f_a,f_b)=H[a,b]
(H[a,a]=0).
"""
import itertools
import numpy as np

INF = 1e300


def nc_matching(pts, dist):
    k = len(pts)
    if k % 2:
        return INF
    if k == 0:
        return 0.0
    M = [[0.0] * (k + 1) for _ in range(k + 1)]
    for length in range(2, k + 1, 2):
        for a in range(0, k - length + 1):
            b = a + length
            best = INF
            pa = pts[a]
            for cc in range(a + 1, b, 2):
                val = dist(pa, pts[cc]) + M[a + 1][cc] + M[cc + 1][b]
                if val < best:
                    best = val
            M[a][b] = best
    return M[0][k]


class Circle:
    def __init__(self, seq, L, C, H, k=4):
        """seq: cyclic list of items: ('arc', label) or ('free', a) or ('edge', i).
        Must alternate sensibly: edges separate labelled arcs; a ('free', a)
        item sits between two labelled arcs and is expanded into k copies."""
        self.L, self.C, self.H = L, C, H
        pts, labels = [], []
        # expand: walk items; each 'arc' item contributes an arc label, 'edge' and 'free' contribute points
        for it in seq:
            if it[0] == 'arc':
                labels.append(it[1])
            elif it[0] == 'edge':
                pts.append(('e', it[1]))
            elif it[0] == 'free':
                for c in range(k):
                    pts.append(('f', it[1]))
                    if c < k - 1:
                        labels.append(None)
        # after the walk, arcs and points must alternate: arc0 pt0 arc1 pt1 ... -> we stored
        # labels in order; arrange so that arc j lies between pts[j-1] and pts[j]
        assert len(pts) == len(labels), (len(pts), len(labels))
        self.pts = pts; self.labels = labels  # arc j precedes point j (cyclic)
        self.free_idx = [j for j, l in enumerate(labels) if l is None]

    def dist(self, p, q):
        if p[0] == 'e' and q[0] == 'e':
            return self.L[p[1], q[1]]
        if p[0] == 'f' and q[0] == 'f':
            return 0.0 if p[1] == q[1] else self.H[p[1], q[1]]
        if p[0] == 'f':
            p, q = q, p
        return self.C[q[1]][p[1]]

    def S(self, A):
        """A: set of party labels (0..5).  min over free memberships."""
        n = len(self.pts)
        base = [(l in A) if l is not None else None for l in self.labels]
        best = INF
        for bits in itertools.product((False, True), repeat=len(self.free_idx)):
            mem = list(base)
            for j, b in zip(self.free_idx, bits):
                mem[j] = b
            # point j between arc j and arc j+1
            act = [self.pts[j] for j in range(n) if mem[j] != mem[(j + 1) % n]]
            v = nc_matching(act, self.dist)
            if v < best:
                best = v
        return best
