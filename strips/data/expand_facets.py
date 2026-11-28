"""Expand the eight exact n=5 facet representatives through purification symmetry.

Input is the author-maintained HEC database's data/n=5/facets.json, saved as
facets_n5_upstream.json in this directory. No third-party dependencies.
"""

from itertools import combinations, permutations
import json
from pathlib import Path


def expand(rows, n):
    subsets = [c for k in range(1, n + 1) for c in combinations(range(n), k)]
    index = {frozenset(c): i for i, c in enumerate(subsets)}
    full = frozenset(range(n + 1))
    orbits = []
    for row in rows:
        assert len(row) == 2**n - 1
        images = set()
        for permutation in permutations(range(n + 1)):
            result = [0] * len(subsets)
            for subset, coefficient in zip(subsets, row):
                if coefficient == 0:
                    continue
                image = frozenset(permutation[j] for j in subset)
                if n in image:
                    image = full - image
                result[index[image]] += coefficient
            images.add(tuple(result))
        orbits.append(sorted(images))
    allrows = sorted(set(row for orbit in orbits for row in orbit))
    lookup = {row: i for i, row in enumerate(allrows)}
    return subsets, orbits, allrows, lookup


def main():
    root = Path(__file__).resolve().parent
    rows = json.loads((root / "facets_n5_upstream.json").read_text())
    subsets, orbits, allrows, lookup = expand(rows, 5)
    assert len(allrows) == 372
    assert sorted(map(len, orbits)) == [10, 15, 20, 45, 60, 60, 72, 90]
    result = {
        "source": "https://github.com/SergioHC95/Holographic-Entropy-Cone/blob/master/data/n=5/facets.json",
        "paper": "https://arxiv.org/abs/1903.09148",
        "n": 5,
        "coefficient_convention": "dot(coefficients,entropy_vector)>=0",
        "subset_order": ["".join("ABCDE"[i] for i in c) for c in subsets],
        "orbit_sizes": [len(orbit) for orbit in orbits],
        "facets": allrows,
        "facet_orbits": [[lookup[row] for row in orbit] for orbit in orbits],
    }
    output = root / "facets_n5_all372.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Wrote {len(allrows)} exact facets to {output}")
    print("Orbit sizes:", result["orbit_sizes"])


if __name__ == "__main__":
    main()
