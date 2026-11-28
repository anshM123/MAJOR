#!/usr/bin/env python3
"""Exact finite-data check of the chord-arrangement mincut theorem.

No third-party packages are required.  Geometry uses Fraction coordinates,
capacities and maxflows are integers.  The check constructs the full arrangement
dual, not the complete Bell-pair graph, and compares mincuts over unrestricted
interior cell colorings against independent noncrossing matching dynamic
programming for every boundary selection (up to complementation).
"""

from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction as Q
from functools import cmp_to_key, lru_cache
import json
from pathlib import Path
import random
import time


def sub(p, q):
    return p[0] - q[0], p[1] - q[1]


def cross(u, v):
    return u[0] * v[1] - u[1] * v[0]


def orient(a, b, c):
    return cross(sub(b, a), sub(c, a))


def between(a, b, p):
    return orient(a, b, p) == 0 and all(
        min(a[k], b[k]) <= p[k] <= max(a[k], b[k]) for k in (0, 1)
    )


def proper_cross(a, b, c, d):
    return orient(a, b, c) * orient(a, b, d) < 0 and orient(c, d, a) * orient(c, d, b) < 0


def intersection(a, b, c, d):
    u, v = sub(b, a), sub(d, c)
    t = cross(sub(c, a), v) / cross(u, v)
    return a[0] + t * u[0], a[1] + t * u[1]


def interpolate(a, b, t):
    return a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])


def edge(a, b):
    return (a, b) if a < b else (b, a)


def split_matrix(w):
    n = len(w)
    D = [[0] * n for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            D[a][b] = D[b][a] = sum(
                w[i][j]
                for i in range(a, b)
                for j in list(range(b, n)) + list(range(a))
            )
    for i in range(n):
        for j in range(i + 1, n):
            reconstructed_twice = (
                D[i][j] + D[(i + 1) % n][(j + 1) % n]
                - D[i][(j + 1) % n] - D[(i + 1) % n][j]
            )
            assert reconstructed_twice == 2 * w[i][j]
    return D


def arrangement(n):
    """Use a rational convex polygon; each exterior cap is one terminal."""
    points = [(Q(i), Q(i * i)) for i in range(n)]
    chords = [(points[i], points[j]) for i in range(n) for j in range(i + 1, n)]
    vertices_on_chord = [{a, b} for a, b in chords]
    for k, (a, b) in enumerate(chords):
        for ell in range(k + 1, len(chords)):
            c, d = chords[ell]
            if proper_cross(a, b, c, d):
                p = intersection(a, b, c, d)
                vertices_on_chord[k].add(p)
                vertices_on_chord[ell].add(p)
    edges = set()
    vertices = set()
    for (a, b), vv in zip(chords, vertices_on_chord):
        vv = sorted(vv, key=lambda p: (p[0] - a[0]) / (b[0] - a[0]))
        vertices.update(vv)
        edges.update(edge(c, d) for c, d in zip(vv, vv[1:]))
    edges = sorted(edges)

    neighbors = {p: [] for p in vertices}
    for a, b in edges:
        neighbors[a].append(b)
        neighbors[b].append(a)

    def angle_comparator(center):
        def half(v):
            return 0 if v[1] > 0 or (v[1] == 0 and v[0] >= 0) else 1

        def compare(p, q):
            u, v = sub(p, center), sub(q, center)
            if half(u) != half(v):
                return -1 if half(u) < half(v) else 1
            c = cross(u, v)
            return -1 if c > 0 else (1 if c < 0 else 0)

        return compare

    for p in vertices:
        neighbors[p].sort(key=cmp_to_key(angle_comparator(p)))

    face_of = {}
    face_polygons = []
    for a, b in edges:
        for start in ((a, b), (b, a)):
            if start in face_of:
                continue
            polygon = []
            h = start
            while h not in face_of:
                face_of[h] = len(face_polygons)
                u, v = h
                polygon.append(u)
                vv = neighbors[v]
                h = v, vv[(vv.index(u) - 1) % len(vv)]
            assert h == start
            face_polygons.append(polygon)

    areas = [sum(cross(a, b) for a, b in zip(pp, pp[1:] + pp[:1])) for pp in face_polygons]
    outer = [i for i, a in enumerate(areas) if a < 0]
    assert len(outer) == 1 and all(a != 0 for a in areas)
    outer = outer[0]
    face_ids = {i: k for k, i in enumerate(i for i in range(len(areas)) if i != outer)}
    boundary = {edge(points[i], points[(i + 1) % n]): i for i in range(n)}

    tests = []
    for i in range(n):
        for j in range(i + 1, n):
            attempt = 0
            while True:
                attempt += 1
                # Fix the first endpoint and vary only the second.  Varying
                # both by a projectively related formula can accidentally
                # keep every trial on a pencil through one forbidden vertex.
                a = interpolate(points[i], points[(i + 1) % n], Q(1, 3))
                b = interpolate(points[j], points[(j + 1) % n], Q(1, attempt + 3))
                if not any(between(a, b, p) for p in vertices):
                    break
                assert attempt < 10000
            tests.append((i, j, a, b))

    # Store the list of test-chord labels crossing each primal edge so one
    # exact geometry can be reused for many capacity assignments.
    dual_edges = []
    n_internal = len(face_ids)
    for a, b in edges:
        left, right = face_of[a, b], face_of[b, a]
        if edge(a, b) in boundary:
            i = boundary[edge(a, b)]
            assert (left == outer) ^ (right == outer)
            inside = right if left == outer else left
            crossed = [(ii, jj) for ii, jj, _, _ in tests if i in (ii, jj)]
            u, v = face_ids[inside], n_internal + i
        else:
            assert left != outer and right != outer
            crossed = [(i, j) for i, j, c, d in tests if proper_cross(a, b, c, d)]
            u, v = face_ids[left], face_ids[right]
        dual_edges.append((u, v, crossed))

    return {
        "n": n,
        "nodes": n_internal + n,
        "internal_nodes": n_internal,
        "terminals": list(range(n_internal, n_internal + n)),
        "dual_edges": dual_edges,
        "arrangement_vertices": len(vertices),
        "arrangement_edges": len(edges),
    }


class Dinic:
    def __init__(self, n):
        self.g = [[] for _ in range(n)]

    def add_edge(self, u, v, c):
        self.g[u].append([v, len(self.g[v]), c])
        self.g[v].append([u, len(self.g[u]) - 1, 0])

    def maxflow(self, source, sink):
        flow = 0
        while True:
            level = [-1] * len(self.g)
            level[source] = 0
            queue = deque([source])
            while queue:
                u = queue.popleft()
                for v, _, c in self.g[u]:
                    if c and level[v] < 0:
                        level[v] = level[u] + 1
                        queue.append(v)
            if level[sink] < 0:
                return flow
            pos = [0] * len(self.g)

            def push(u, amount):
                if u == sink:
                    return amount
                while pos[u] < len(self.g[u]):
                    e = self.g[u][pos[u]]
                    v, reverse, capacity = e
                    if capacity and level[v] == level[u] + 1:
                        sent = push(v, min(amount, capacity))
                        if sent:
                            e[2] -= sent
                            self.g[v][reverse][2] += sent
                            return sent
                    pos[u] += 1
                return 0

            while True:
                sent = push(source, sum(c for _, _, c in self.g[source]) + 1)
                if not sent:
                    break
                flow += sent


def weighted_edges(geo, w):
    return [(u, v, sum(w[i][j] for i, j in crossed)) for u, v, crossed in geo["dual_edges"]]


def cut_value(geo, weighted, mask):
    source, sink = geo["nodes"], geo["nodes"] + 1
    flow = Dinic(geo["nodes"] + 2)
    high = 1 + sum(c for _, _, c in weighted)
    for u, v, c in weighted:
        flow.add_edge(u, v, c)
        flow.add_edge(v, u, c)
    for i, terminal in enumerate(geo["terminals"]):
        if mask & (1 << i):
            flow.add_edge(source, terminal, high)
        else:
            flow.add_edge(terminal, sink, high)
    return flow.maxflow(source, sink)


def matching_value(D, mask):
    n = len(D)
    endpoints = tuple(i for i in range(n) if bool(mask & (1 << i)) != bool(mask & (1 << ((i - 1) % n))))

    @lru_cache(None)
    def solve(xx):
        if not xx:
            return 0
        return min(
            D[xx[0]][xx[j]] + solve(xx[1:j]) + solve(xx[j + 1 :])
            for j in range(1, len(xx), 2)
        )

    return solve(endpoints)


def brute_cut_value(geo, weighted, mask):
    """Independent exhaustive check of every interior cell coloring."""
    best = None
    for inner in range(1 << geo["internal_nodes"]):
        selected = inner | (mask << geo["internal_nodes"])
        value = sum(c for u, v, c in weighted if bool(selected & (1 << u)) != bool(selected & (1 << v)))
        best = value if best is None else min(best, value)
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert 3 <= args.max_n <= 12
    rng = random.Random(20261128)
    report = {"arithmetic": "exact rational geometry; integer capacities and flows", "cases": []}
    total = 0
    start = time.monotonic()
    for n in range(3, args.max_n + 1):
        geo = arrangement(n)
        cases = [
            ("all_one", [[int(i != j) for j in range(n)] for i in range(n)]),
            ("random_with_zero", [[0] * n for _ in range(n)]),
            ("sparse_crossing", [[0] * n for _ in range(n)]),
        ]
        for i in range(n):
            for j in range(i + 1, n):
                cases[1][1][i][j] = cases[1][1][j][i] = rng.randrange(0, 10)
        cases[2][1][0][n // 2] = cases[2][1][n // 2][0] = 5
        cases[2][1][1][(n // 2 + 1) % n] = cases[2][1][(n // 2 + 1) % n][1] = 7
        for name, w in cases:
            D = split_matrix(w)
            weighted = weighted_edges(geo, w)
            checked = 0
            for mask in range(1 << (n - 1)):
                expected = matching_value(D, mask)
                observed = cut_value(geo, weighted, mask)
                assert observed == expected, (n, name, mask, expected, observed)
                if n <= 4:
                    exhaustive = brute_cut_value(geo, weighted, mask)
                    assert observed == exhaustive
                checked += 1
            total += checked
            report["cases"].append({
                "endpoints": n, "weights": name, "graph_nodes": geo["nodes"],
                "graph_edges": len(weighted), "boundary_assignments": checked,
                "exhaustive_interior_colorings": n <= 4,
                "status": "all exact equalities passed",
            })
        print(f"n={n}: {len(cases) * (1 << (n - 1))} exact mincut/matching equalities; {geo['nodes']} graph nodes", flush=True)

    # The full Bell-pair graph fails for a disconnected selection; the
    # arrangement graph captures the missing minimum over internal topology.
    w = [[0] * 4 for _ in range(4)]
    w[0][2] = w[2][0] = 5
    w[1][3] = w[3][1] = 7
    D = split_matrix(w)
    mask = 0b0101
    bell = sum(w[i][j] for i in range(4) for j in range(i + 1, 4) if bool(mask & (1 << i)) != bool(mask & (1 << j)))
    arrangement_result = cut_value(arrangement(4), weighted_edges(arrangement(4), w), mask)
    expected = matching_value(D, mask)
    assert bell == 0 and arrangement_result == expected == 10
    report["bell_pair_countercheck"] = {"bell_pair_graph_cut": bell, "arrangement_graph_cut": arrangement_result, "matching_minimum": expected}
    report["total_exact_equalities"] = total
    report["seconds"] = round(time.monotonic() - start, 3)
    report["status"] = "PASS"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "cases"}, indent=2))


if __name__ == "__main__":
    main()
