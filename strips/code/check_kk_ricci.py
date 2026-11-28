"""Exact local checks of the KK horizontal Ricci identity.

This is a sign/coefficient check, not the proof of the identity.  It evaluates
the full coordinate Ricci tensor from metric two-jets, without using a KK
reduction identity, and compares its horizontal block with the stated formula.
"""

import sympy as s

t, x, z = base = s.symbols("t x z", real=True)
eta = s.diag(-1, 1, 1)
origin = dict.fromkeys(base, 0)


def at_origin(expr):
    return s.simplify(expr.subs(origin))


def ricci_at_origin(metric, coordinates):
    dim = len(coordinates)
    g = metric.applyfunc(at_origin)
    gi = g.inv()
    dg = [metric.diff(q).applyfunc(at_origin) for q in coordinates]
    ddg = [[metric.diff(q, r).applyfunc(at_origin)
            for r in coordinates] for q in coordinates]
    dgi = [-gi * dg[q] * gi for q in range(dim)]
    gamma = s.MutableDenseNDimArray.zeros(dim, dim, dim)
    dgamma = s.MutableDenseNDimArray.zeros(dim, dim, dim, dim)
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                gamma[a, b, c] = sum(
                    gi[a, r] * (dg[b][r, c] + dg[c][r, b] - dg[r][b, c])
                    for r in range(dim)) / 2
                for q in range(dim):
                    dgamma[q, a, b, c] = sum(
                        dgi[q][a, r] *
                        (dg[b][r, c] + dg[c][r, b] - dg[r][b, c])
                        + gi[a, r] * (ddg[q][b][r, c]
                                     + ddg[q][c][r, b]
                                     - ddg[q][r][b, c])
                        for r in range(dim)) / 2
    result = s.zeros(dim)
    for a in range(dim):
        for b in range(dim):
            result[a, b] = s.simplify(sum(
                dgamma[c, c, a, b] - dgamma[b, c, a, c]
                + sum(gamma[c, c, r] * gamma[r, a, b]
                      - gamma[c, b, r] * gamma[r, a, c]
                      for r in range(dim))
                for c in range(dim)))
    return result


def check(label, G, A):
    m = G.rows
    ys = s.symbols(f"y0:{m}")
    metric = (eta + A.T * G * A).row_join(A.T * G)
    metric = metric.col_join((G * A).row_join(G))
    ricci = ricci_at_origin(metric, (*base, *ys))
    horizontal = s.eye(3).col_join(-A.applyfunc(at_origin))
    rh = horizontal.T * ricci * horizontal
    logvolume = s.log(G.det()) / 2
    G0 = G.applyfunc(at_origin)
    P = [(G.inv() * G.diff(q)).applyfunc(at_origin) for q in base]
    F = [s.Matrix(3, 3, lambda a, b:
                  at_origin(s.diff(A[i, b], base[a])
                            - s.diff(A[i, a], base[b])))
         for i in range(m)]
    predicted = s.zeros(3)
    for a in range(3):
        for b in range(3):
            predicted[a, b] = (
                -at_origin(s.diff(logvolume, base[a], base[b]))
                -s.trace(P[a] * P[b]) / 4
                -sum(G0[i, j] * sum(F[i][a, c] * eta[c, r] * F[j][b, r]
                                   for c in range(3) for r in range(3))
                     for i in range(m) for j in range(m)) / 2)
    residual = s.simplify(rh - predicted)
    assert residual == s.zeros(3), (label, residual)
    ge = G.det() * eta
    re = ricci_at_origin(ge, base)
    k = s.Matrix([1, s.Rational(3, 5), s.Rational(4, 5)])
    ku = sum(k[a] * at_origin(s.diff(logvolume, base[a])) for a in range(3))
    pk = sum((k[a] * P[a] for a in range(3)), s.zeros(m))
    fterm = sum(G0[i, j] * (k.T * F[i] * eta * F[j].T * k)[0]
                for i in range(m) for j in range(m)) / 2
    positive = ku ** 2 + s.trace(pk * pk) / 4 + fterm
    null_residual = s.simplify((k.T * (re - rh) * k)[0] - positive)
    assert null_residual == 0, (label, null_residual)
    print(f"{label}: horizontal Ricci residual is zero; null residual is zero; "
          f"positive correction = {positive}")


check("m=1 with radius and twist",
      s.Matrix([[s.exp(2 * (t + 2*x + 3*z + t*x))]]),
      s.Matrix([[6*z + 1, 4*t + 5*z + 2, t*x + 3]]))

L = s.Matrix([[1, t + 2*x - z], [0, 1]])
G2 = L.T * s.diag(s.exp(2 * (t + x + t*z)),
                  s.exp(2 * (x + 2*z + x*z))) * L
A2 = s.Matrix([[z + 1, t + 2*z + 2, x + 3],
               [2*x + 4, 3*t + 5, t - x + 6]])
check("m=2 with shape, volume, and two twists", G2, A2)
