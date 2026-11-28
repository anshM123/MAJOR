#!/usr/bin/env python3
"""Reproduce symbolic obstruction and finite-cut BTZ limit checks.

The numerical checks illustrate, and do not replace, the estimates proved in
CUSPS_EXTREMAL_BTZ.md.  Only NumPy and SymPy are required.
"""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp


def symbolic_checks():
    lam = sp.symbols("lam", positive=True)
    shift = sp.symbols("shift", nonzero=True)
    entries = sp.symbols("x11 x12 x21 x22")
    X = sp.Matrix(2, 2, entries)
    L = sp.diag(lam, 1 / lam)
    R = sp.Matrix([[1, shift], [0, 1]])
    TX = L * X * R.inv()
    action = sp.Matrix(list(TX)).jacobian(entries)
    a, b = sp.symbols("a b")
    invariant_planes = [
        sp.Matrix([[a, b], [0, 0]]),
        sp.Matrix([[0, 0], [a, b]]),
        sp.Matrix([[0, a], [0, b]]),
    ]
    plane_determinants = [sp.factor(W.det()) for W in invariant_planes]
    assert plane_determinants == [0, 0, 0]

    r, r0 = sp.symbols("r r0", positive=True)
    # Coordinates v, r, psi; only the radial derivative is nonzero.
    metric = sp.Matrix([[-r*r + 2*r0*r0, 1, -r0*r0],
                        [1, 0, 0], [-r0*r0, 0, r*r]])
    inverse = sp.simplify(metric.inv())
    acceleration = [sp.simplify(-r * inverse[mu, 1]).subs(r, r0)
                    for mu in range(3)]
    assert acceleration == [-r0, 0, -r0]
    return {
        "holonomy_action": str(action),
        "invariant_plane_determinants": list(map(str, plane_determinants)),
        "horizon_acceleration_v_r_psi": list(map(str, acceleration)),
        "symbolic_status": "exact checks passed",
    }


def matchings(indices):
    if not indices:
        yield ()
        return
    i = indices[0]
    for j in indices[1:]:
        remaining = tuple(x for x in indices if x != i and x != j)
        for tail in matchings(remaining):
            yield ((i, j),) + tail


def interlaces(a, b, c, d):
    return a < c < b < d or c < a < d < b


def compatible(arcs, period):
    for (a, b), (c, d) in itertools.combinations(arcs, 2):
        for m in range(-2, 3):
            if interlaces(a, b, c + m*period, d + m*period):
                return False
    return True


def candidates(points, period):
    """All embedded reduced cuts, outer trace zero at phi=0.

    An arc stores (endpoint index, endpoint index, winding), with winding
    zero or minus one.  The core is included iff an odd number of its
    disk-side boundary intervals wrap through phi=0.
    """
    for pairs in matchings(tuple(range(len(points)))):
        for wraps in itertools.product((0, 1), repeat=len(pairs)):
            intervals = []
            components = []
            for (i, j), wrap in zip(pairs, wraps):
                assert i < j
                p, q = points[i], points[j]
                intervals.append((q, p + period) if wrap else (p, q))
                components.append((i, j, -wrap))
            if compatible(intervals, period):
                yield tuple(components), sum(wraps) % 2


def log_sinh_over_parameter(parameter, width):
    if parameter == 0:
        return math.log(width)
    x = parameter * width
    return math.log(math.sinh(x) / parameter)


def check_btz_family():
    period = 2*math.pi
    beta = 0.73
    epsilon = 1e-4
    parties = [(0.20, 0.90), (1.00, 1.80), (1.95, 2.75),
               (2.90, 3.70), (3.85, 4.90)]
    alphas = [0.0, 0.003, 0.01, 0.03, 0.1, 0.3]
    # The derivative is bounded by 0.30, giving a periodic spacelike curve.
    boundary_time = lambda x: 0.20*math.sin(x) + 0.05*math.cos(2*x)
    entropies = {alpha: {0: 0.0} for alpha in alphas}
    counts = {}
    error_rows = []
    for mask in range(1, 1 << len(parties)):
        points = sorted(x for p, interval in enumerate(parties)
                        if mask & (1 << p) for x in interval)
        times = [boundary_time(x) for x in points]
        cuts = list(candidates(points, period))
        k = len(points)//2
        counts.setdefault(str(k), []).append(len(cuts))
        for alpha in alphas:
            def cut_cost(cut):
                components, has_core = cut
                value = has_core * period * (alpha + beta)
                for i, j, winding in components:
                    u = abs(points[j] - points[i] + times[j] - times[i]
                            + winding*period)
                    v = abs(points[j] - points[i] - times[j] + times[i]
                            + winding*period)
                    assert 0 < u < period and 0 < v < period
                    value += (log_sinh_over_parameter(alpha, u)
                              + log_sinh_over_parameter(beta, v)
                              - 2*math.log(epsilon))
                return value
            entropies[alpha][mask] = min(map(cut_cost, cuts))
        for alpha in alphas[1:]:
            difference = entropies[alpha][mask] - entropies[0.0][mask]
            bound = k*alpha**2*period**2/6 + period*alpha
            assert -1e-12 <= difference <= bound + 1e-12
            error_rows.append({"mask": mask, "alpha": alpha,
                               "difference": difference, "bound": bound})

    # Independent sample checks of two important balanced cone inequalities.
    facet_rows = []
    for alpha, values in entropies.items():
        mmi = []
        for a, b, c in itertools.combinations(range(1, 32), 3):
            if a & b or a & c or b & c:
                continue
            mmi.append(values[a | b] + values[a | c] + values[b | c]
                       - values[a] - values[b] - values[c] - values[a | b | c])
        cyclic = (sum(values[sum(1 << ((i+j) % 5) for j in range(3))]
                      for i in range(5))
                  - sum(values[(1 << i) | (1 << ((i+1) % 5))]
                        for i in range(5)) - values[31])
        assert min(mmi) >= -1e-10 and cyclic >= -1e-10
        facet_rows.append({"alpha": alpha, "minimum_MMI": min(mmi),
                           "five_party_cyclic": cyclic})

    # Exhibit a genuine limiting core branch for one large interval.
    ell = 0.99*period
    direct = (math.log(ell) + math.log(math.sinh(beta*ell)/beta)
              - 2*math.log(epsilon))
    complement_core = (math.log(period-ell)
                       + math.log(math.sinh(beta*(period-ell))/beta)
                       - 2*math.log(epsilon) + period*beta)
    assert complement_core < direct
    return {
        "normalization": "AdS radius=1; length, before division by 4G_N",
        "period": period, "beta": beta, "epsilon": epsilon,
        "party_intervals": parties,
        "time_profile": "0.20*sin(phi)+0.05*cos(2*phi)",
        "candidate_counts_by_number_of_intervals":
            {k: sorted(set(v)) for k, v in counts.items()},
        "uniform_estimate_checks": len(error_rows),
        "largest_observed_error_to_bound_ratio":
            max(row["difference"]/row["bound"] for row in error_rows),
        "selected_inequality_checks": facet_rows,
        "large_interval_core_branch": {
            "width": ell, "direct_length": direct,
            "complement_plus_limiting_core_length": complement_core,
            "limiting_core_is_an_attained_geodesic": False,
        },
        "entropy_vectors": {str(a): v for a, v in entropies.items()},
        "numerical_status": "all stated checks passed",
    }


if __name__ == "__main__":
    output = symbolic_checks() | check_btz_family()
    destination = Path(__file__).with_name("cusp_validation.json")
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: value for key, value in output.items()
                      if key != "entropy_vectors"}, indent=2))
