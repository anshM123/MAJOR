#!/usr/bin/env python3
"""Exact finite-group mincut/quotient check; integer cover, rational quotient.

This verifies the combinatorial descent lemma, not the spacetime hypotheses.
Only the Python standard library is required.
"""

from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import random


def orbit(seed, perm):
    out = []
    point = seed
    while point not in out:
        out.append(point)
        point = perm[point]
    return tuple(out)


def edge_orbit(edge, perm):
    out = set()
    point = edge
    while point not in out:
        out.add(point)
        point = tuple(sorted((perm[point[0]], perm[point[1]])))
    return tuple(sorted(out))


def vertex_orbits(perm):
    seen = set()
    out = []
    for point in range(len(perm)):
        if point not in seen:
            cyc = orbit(point, perm)
            out.append(cyc)
            seen.update(cyc)
    return out


def capacity(selected, edges):
    return sum(w for a, b, w in edges if (a in selected) != (b in selected))


def mincuts(vertices, edges, ones, zeros):
    free = sorted(set(vertices) - ones - zeros)
    best = None
    minimizers = []
    for bits in range(1 << len(free)):
        chosen = ones | {v for k, v in enumerate(free) if bits >> k & 1}
        val = capacity(chosen, edges)
        if best is None or val < best:
            best, minimizers = val, [chosen]
        elif val == best:
            minimizers.append(chosen)
    return best, minimizers


def verify_case(q, perm, edges, ones, zeros):
    cycles = vertex_orbits(perm)
    assert all(q % len(cyc) == 0 for cyc in cycles)
    label = {v: k for k, cyc in enumerate(cycles) for v in cyc}
    assert {perm[v] for v in ones} == ones
    assert {perm[v] for v in zeros} == zeros
    weights = {(a, b): w for a, b, w in edges}
    for a, b, w in edges:
        assert weights[tuple(sorted((perm[a], perm[b])))] == w

    best, minimizers = mincuts(range(len(perm)), edges, ones, zeros)
    noninvariant = 0
    for chosen in minimizers:
        if {perm[v] for v in chosen} != chosen:
            noninvariant += 1
        closure = set(chosen)
        for _ in range(q):
            closure.update(perm[v] for v in tuple(closure))
        assert capacity(closure, edges) == best
        assert {perm[v] for v in closure} == closure

    quotient_weights = {}
    for a, b, w in edges:
        i, j = sorted((label[a], label[b]))
        if i != j:
            quotient_weights[i, j] = quotient_weights.get((i, j), 0) + Fraction(w, q)
    quotient_edges = [(i, j, w) for (i, j), w in quotient_weights.items()]
    quotient_ones = {label[v] for v in ones}
    quotient_zeros = {label[v] for v in zeros}
    quotient_best, _ = mincuts(range(len(cycles)), quotient_edges, quotient_ones, quotient_zeros)
    assert Fraction(best, q) == quotient_best
    invariant_best = min(
        capacity(set().union(*(set(cycles[k]) for k in selected)), edges)
        for bits in range(1 << len(cycles))
        if (selected := {k for k in range(len(cycles)) if bits >> k & 1}) >= quotient_ones
        and not selected & quotient_zeros
    )
    assert invariant_best == best
    return {
        "vertices": len(perm),
        "group_order": q,
        "cover_minimum": str(best),
        "quotient_minimum": str(quotient_best),
        "minimum_count": len(minimizers),
        "noninvariant_minimum_count": noninvariant,
        "has_fixed_vertices": any(len(cyc) == 1 for cyc in cycles),
        "has_stabilized_edges": any(len(edge_orbit((a, b), perm)) < q for a, b, _ in edges),
    }


def random_case(q, trial, rng):
    # Three fixed terminal vertices plus mixed-length internal cycles.
    divisors = [k for k in range(1, q + 1) if q % k == 0]
    sizes = [1, 1, 1, q]
    spare = 12 - sum(sizes)
    while spare > 0:
        options = [k for k in divisors if k <= spare]
        k = rng.choice(options)
        sizes.append(k)
        spare -= k
        if rng.random() < 0.35:
            break
    perm = []
    for size in sizes:
        start = len(perm)
        perm.extend(start + (j + 1) % size for j in range(size))
    weights = {}
    for edge in combinations(range(len(perm)), 2):
        if edge not in weights:
            w = rng.randrange(0, 6)
            for translated in edge_orbit(edge, perm):
                weights[translated] = w
    edges = [(a, b, w) for (a, b), w in weights.items()]
    mask = 1 + trial % 6
    ones = {v for v in range(3) if mask >> v & 1}
    zeros = set(range(3)) - ones
    return verify_case(q, perm, edges, ones, zeros)


def main():
    rng = random.Random(20261128)
    cases = []
    # Explicit degeneracy: the two swapped internal vertices may be assigned
    # independently, giving asymmetric minima whose union is symmetric.
    cases.append(verify_case(
        2, [0, 1, 3, 2],
        [(0, 1, 0), (0, 2, 1), (0, 3, 1), (1, 2, 1), (1, 3, 1), (2, 3, 0)],
        {0}, {1},
    ))
    for q in (2, 3, 4, 6):
        cases.extend(random_case(q, trial, rng) for trial in range(80))
    result = {
        "status": "all exact checks passed",
        "seed": 20261128,
        "cases": len(cases),
        "group_orders": sorted({case["group_order"] for case in cases}),
        "cases_with_noninvariant_minima": sum(case["noninvariant_minimum_count"] > 0 for case in cases),
        "cases_with_stabilized_edges": sum(case["has_stabilized_edges"] for case in cases),
        "explicit_degenerate_case": cases[0],
        "interpretation": "Checks submodular invariant-minimum descent and capacity normalization only; not geometric or physical hypotheses.",
    }
    destination = Path(__file__).with_suffix(".json")
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
