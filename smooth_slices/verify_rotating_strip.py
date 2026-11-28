"""Exact rotating-strip algebra for SMOOTH_SLICES.md §9; no sampled inequalities."""

import json
from pathlib import Path

import sympy as S

s, r, sp, rp, q, tau = S.symbols("s r sp rp q tau", real=True)
a, b, z = S.symbols("a b z", positive=True)


def N(X, Y):
    return X[0, 0] * Y[1, 1] + X[1, 1] * Y[0, 0] - X[0, 1] * Y[1, 0] - X[1, 0] * Y[0, 1]


def metric(X, Y):
    return -N(X, Y) / 2


def zero(expr):
    if S.simplify(S.powsimp(expr, force=True)) == 0:
        return True
    return S.simplify(S.expand(expr.rewrite(S.exp))) == 0


checks = {}


def check(name, expr):
    expressions = list(expr) if isinstance(expr, S.MatrixBase) else [expr]
    checks[name] = all(zero(item) for item in expressions)
    if not checks[name]:
        raise AssertionError(f"Exact check failed: {name}: {expr}")


c = S.Matrix([[0, S.exp(s)], [-S.exp(-s), 0]])
n = S.diag(S.exp(q * s), -S.exp(-q * s))
T = S.diag(S.exp(q * s), S.exp(-q * s))
e = c.diff(s)
frame = [c, e, n, T]
gram = S.Matrix([[metric(u, v) for v in frame] for u in frame])
check("ambient_frame_gram", gram - S.diag(-1, 1, 1, -1))
check("frame_derivatives", S.Matrix.hstack(c.diff(s) - e, e.diff(s) - c, n.diff(s) - q * T, T.diff(s) - q * n))

Phi = S.cosh(r) * c + S.sinh(r) * n
check("Phi_determinant", Phi.det() - 1)
Pr, Ps = Phi.diff(r), Phi.diff(s)
check("Phi_metric_rr", metric(Pr, Pr) - 1)
check("Phi_metric_rs", metric(Pr, Ps))
check("Phi_metric_ss", metric(Ps, Ps) - (1 + (1 - q**2) * S.sinh(r)**2))
check("T_tangent_to_AdS", metric(T, Phi))

Phi2 = Phi.subs({s: sp, r: rp}, simultaneous=True)
delta = s - sp
pair = S.cosh(r) * S.cosh(rp) * S.cosh(delta) - S.sinh(r) * S.sinh(rp) * S.cosh(q * delta)
check("pairwise_N_formula", N(Phi, Phi2) / 2 - pair)
positive_form = S.cosh(r - rp) * S.cosh(q * delta) + S.cosh(r) * S.cosh(rp) * (S.cosh(delta) - S.cosh(q * delta))
check("pairwise_positive_decomposition", pair - positive_form)
check("positive_conformal_sine", (Phi[0, 1] - Phi[1, 0]) / 2 - S.cosh(r) * S.cosh(s))

L = (a + b) / 2
q_hol = (a - b) / (a + b)
hL = S.diag(S.exp(a / 2), S.exp(-a / 2))
hR_inverse = S.diag(S.exp(-b / 2), S.exp(b / 2))
Phi_hol = Phi.subs(q, q_hol)
check("holonomy_equivariance", hL * Phi_hol * hR_inverse - Phi_hol.subs(s, s + L))
check("uniform_margin_holonomy_formula", (1 - q_hol**2) / (1 + q_hol**2) - 2 * a * b / (a**2 + b**2))

alpha, beta = S.symbols("alpha beta", real=True)
v = alpha * Pr + beta * Ps
check("auxiliary_Riemannian_metric", metric(v, v) + 2 * metric(T, v)**2 - alpha**2 - (S.cosh(r)**2 + q**2 * S.sinh(r)**2) * beta**2)
normal_numerator = S.cosh(r) * T + q * S.sinh(r) * e
check("normal_orthogonal_r", metric(normal_numerator, Pr))
check("normal_orthogonal_s", metric(normal_numerator, Ps))
check("normal_norm", metric(normal_numerator, normal_numerator) + S.cosh(r)**2 - q**2 * S.sinh(r)**2)

F = S.cosh(r) * (S.cos(tau) * c + S.sin(tau) * T) + S.sinh(r) * n
check("F_determinant", F.det() - 1)
coordinates = [r, s, tau]
derivatives = [F.diff(x) for x in coordinates]
F_metric = S.Matrix([[metric(u, v) for v in derivatives] for u in derivatives])
expected = S.Matrix([
    [1, q * S.sin(tau), 0],
    [q * S.sin(tau), q**2 + (1 - q**2) * S.cosh(r)**2 * S.cos(tau)**2, -q * S.sinh(r) * S.cosh(r) * S.cos(tau)],
    [0, -q * S.sinh(r) * S.cosh(r) * S.cos(tau), -S.cosh(r)**2],
])
check("F_full_metric", F_metric - expected)
check("F_metric_determinant", expected.det() + S.cosh(r)**4 * S.cos(tau)**2)
cofactor_tau_tau = expected[:2, :2].det()
check("inverse_metric_tau_tau", cofactor_tau_tau - S.cos(tau)**2 * (q**2 + (1 - q**2) * S.cosh(r)**2))

cosh_end, sinh_end = 1 / z + z / 4, 1 / z - z / 4
compact_ss = (1 - q**2) + (1 + q**2) * z**2 / 2 + (1 - q**2) * z**4 / 16
check("compactified_induced_metric", z**2 * (cosh_end**2 - q**2 * sinh_end**2) - compact_ss)
# dr=-dz/z at the positive end, followed by multiplication by z^2.
J = S.diag(-1 / z, 1, 1)
compact_g = z**2 * J.T * expected.subs({S.cosh(r): cosh_end, S.sinh(r): sinh_end}) * J
boundary_g = compact_g.applyfunc(lambda x: S.limit(x, z, 0))
expected_boundary = S.Matrix([
    [1, 0, 0],
    [0, (1 - q**2) * S.cos(tau)**2, -q * S.cos(tau)],
    [0, -q * S.cos(tau), -1],
])
check("compactified_ambient_boundary_metric", boundary_g - expected_boundary)
check("compactified_ambient_nondegenerate", expected_boundary.det() + S.cos(tau)**2)

out = {
    "purpose": "Exact local algebra audit; inequalities are proved analytically in SMOOTH_SLICES.md §9",
    "sympy_version": S.__version__,
    "checks": checks,
    "all_passed": all(checks.values()),
    "number_of_checks": len(checks),
}
destination = Path(__file__).with_name("rotating_strip_verification.json")
destination.write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
