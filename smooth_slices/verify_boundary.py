#!/usr/bin/env python3
"""Exact local identities for SMOOTH_SLICES.md §§2, 4, and 10.

Requires SymPy. No sampling, numerical tolerances, or global causal claims.
Run: python verify_boundary.py [--output boundary_verification.json]
"""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import sympy as sp


def verify() -> dict:
    checks: list[dict] = []

    def zero(name: str, expression) -> None:
        entries = list(expression) if isinstance(expression, sp.MatrixBase) else [expression]
        residuals = [sp.simplify(sp.factor(value)) for value in entries]
        passed = all(value == 0 for value in residuals)
        record = {"name": name, "passed": passed}
        if not passed:
            record["nonzero_residuals"] = [str(value) for value in residuals if value != 0]
        checks.append(record)

    # Generic rank-one endpoint matrices; no special frame is imposed.
    a1, a2, b1, b2, c1, c2, d1, d2 = sp.symbols(
        "a1 a2 b1 b2 c1 c2 d1 d2", real=True
    )
    p = sp.Matrix([a1, a2]) * sp.Matrix([[b1, b2]])
    q = sp.Matrix([c1, c2]) * sp.Matrix([[d1, d2]])
    u = sp.symbols("u", real=True)
    N = sp.expand(sp.trace(p.adjugate() * q))
    X = p + u**2 * q
    w = u * sp.sqrt(N)
    zero("ambient.endpoint_p_rank_one", p.det())
    zero("ambient.endpoint_q_rank_one", q.det())
    zero("ambient.polarization", (p + q).det() - N)
    zero("ambient.null_conic", -X.det() + w**2)
    zero("ambient.tangent_norm_equals_N", -X.diff(u).det() + sp.diff(w, u)**2 - N)
    zero("ambient.reflection_matrix_part", X.subs(u, -u) - X)
    zero("ambient.reflection_normal_part", w.subs(u, -u) + w)
    zero("ambient.endpoint_radial_tangent", X.diff(u).subs(u, 0))

    # Exact Lorentz-boosted Poincare semicircle near the first endpoint.
    z, dx, dt = sp.symbols("z X T", real=True)
    distance2 = dx**2 - dt**2
    branch = (1 - sp.sqrt(1 - 4*z**2/distance2))/2
    x_arc, t_arc = dx*branch, dt*branch
    zero("poincare.endpoint_plane", dt*x_arc - dx*t_arc)
    zero(
        "poincare.endpoint_semicircle",
        (x_arc-dx/2)**2 - (t_arc-dt/2)**2 + z**2 - distance2/4,
    )
    zero("poincare.even_position", sp.Matrix([x_arc, t_arc]).subs(z, -z)
         - sp.Matrix([x_arc, t_arc]))
    zero("poincare.endpoint_position", sp.Matrix([x_arc, t_arc]).subs(z, 0))
    zero("poincare.endpoint_radial_tangent",
         sp.Matrix([sp.diff(x_arc, z), sp.diff(t_arc, z)]).subs(z, 0))
    zero(
        "poincare.conformal_tangent_norm",
        1 + sp.diff(x_arc, z)**2 - sp.diff(t_arc, z)**2
        - distance2/(distance2-4*z**2),
    )

    # Transform the physical BTZ metric from r to z=1/r, then multiply by z^2.
    r, rplus, rminus = sp.symbols("r rplus rminus", real=True)
    lapse2 = (r**2-rplus**2)*(r**2-rminus**2)/r**2
    shift = -rplus*rminus/r**2
    physical_r = sp.Matrix([
        [1/lapse2, 0, 0],
        [0, -lapse2+r**2*shift**2, r**2*shift],
        [0, r**2*shift, r**2],
    ])
    jacobian = sp.diag(-1/z**2, 1, 1)
    transformed = z**2 * jacobian.T * physical_r.subs(r, 1/z) * jacobian
    A = 1/((1-rplus**2*z**2)*(1-rminus**2*z**2))
    h = sp.Matrix([
        [-1+(rplus**2+rminus**2)*z**2, -rplus*rminus*z**2],
        [-rplus*rminus*z**2, 1],
    ])
    compactified = sp.diag(A, 1, 1)
    compactified[1:3, 1:3] = h
    zero("btz.exact_conformal_metric", transformed-compactified)
    zero("btz.even_metric_coefficients", compactified.subs(z, -z)-compactified)
    reflection_jacobian = sp.diag(-1, 1, 1)
    zero("btz.reflection_isometry", reflection_jacobian.T
         * compactified.subs(z, -z) * reflection_jacobian - compactified)
    zero("btz.boundary_metric", compactified.subs(z, 0)-sp.diag(1, -1, 1))
    zero("btz.tangential_determinant", h.det()+1/A)

    # First integrals for a physical unit spacelike geodesic.
    pt, pphi = sp.symbols("p_t p_phi", real=True)
    momentum = sp.Matrix([pt, pphi])
    raised_momentum = h.inv()*momentum
    H = sp.factor((momentum.T*raised_momentum)[0])
    denominator = 1-z**2*H
    ydot = z**2*raised_momentum
    zdot_squared = z**2*denominator/A
    yprime = z*sp.sqrt(A)*raised_momentum/sp.sqrt(denominator)
    zero("ode.conserved_momenta", z**(-2)*h*ydot-momentum)
    zero("ode.physical_unit_speed", A*zdot_squared/z**2+z**2*H-1)
    zero("ode.odd_endpoint_derivative", yprime.subs(z, -z)+yprime)
    zero("ode.conformal_tangent_norm",
         A+(yprime.T*h*yprime)[0]-A/denominator)

    return {
        "all_passed": all(item["passed"] for item in checks),
        "check_count": len(checks),
        "method": "exact symbolic identities; every residual simplifies to zero",
        "python_version": platform.python_version(),
        "sympy_version": sp.__version__,
        "checks": checks,
        "positivity_assumptions": [
            "Ambient endpoint polarization N>0.",
            "Poincare endpoints obey X^2-T^2>0 and |z|<sqrt(X^2-T^2)/2.",
            "In a sufficiently thin BTZ collar A>0 and 1-z^2 H>0.",
        ],
        "scope": (
            "These checks verify local identities only. Embeddedness, compactness, "
            "global acausality, Cauchy extension, and homology require the analytic proofs."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=Path(__file__).with_name("boundary_verification.json"),
    )
    args = parser.parse_args()
    result = verify()
    output = json.dumps(result, indent=2) + "\n"
    args.output.write_text(output, encoding="utf-8")
    print(output, end="")
    raise SystemExit(0 if result["all_passed"] else 1)


if __name__ == "__main__":
    main()
