"""P1: what positivity hypotheses does the Einstein reference itself satisfy, and what do they give for E(0) = 4 pi eta/s
from E(1..N) alone?

Data: E_E(n), n = 1..8 from ../cft_etas/out_einstein_ref.txt (units 2 pi T = 1.5, eta = 1, s = 4 pi, so E_E(0) = 1).
E(n) = int_0^inf Phi(t) e^{-y_n t} dt (times a constant), y_n = 1.5 n.

(H) "Phi >= 0"  <=>  E completely monotone in s (Bernstein)  <=>  E(n) = int_[0,1] x^n dnu(x), nu >= 0, x = e^{-1.5 t}
    (Hausdorff).  Then E(0) = nu([0,1]) (no atom at x = 0 since Phi is integrable... an atom at 0 is invisible anyway).
    LP: min nu([0,1]) subject to the N moments  ->  a rigorous LOWER bound on 4 pi eta/s from Euclidean data alone.
    Upper bound needs a bound int x^{-gamma} dnu <= M (analyticity in Re s > -gamma); we report it for M = M_E-type values.
(S) "Phi completely monotone"  <=>  E Stieltjes: E(s) = int dmu(lam)/(lam + s), mu >= 0.  LP feasibility test with
    residual minimisation.  Einstein's lowest QNM is complex, so (S) must fail; we show at which N the data detect it.

Output: out_p1_einstein_positivity.txt
"""
import numpy as np
from scipy.optimize import linprog

E = np.array([0.4127521731, 0.2381977379, 0.1637545371, 0.1239990926, 0.0995749628,
              0.0831246905, 0.0713148684, 0.0624329775])
TOL = 2e-10


def hausdorff_bounds(N, gamma=None, M=None, tol=TOL):
    x = np.unique(np.concatenate([np.linspace(0, 1, 4001), np.logspace(-12, 0, 2000)]))
    n = np.arange(1, N + 1)[:, None]
    A = x[None, :] ** n
    Aub = np.vstack([A, -A])
    bub = np.concatenate([E[:N] + tol, -(E[:N] - tol)])
    if gamma is not None:
        keep = x > 0
        x, A, Aub = x[keep], A[:, keep], Aub[:, keep]
        w = x ** (-gamma)
        Aub = np.vstack([Aub, w[None, :]])
        bub = np.concatenate([bub, [M]])
    lo = linprog(np.ones_like(x), A_ub=Aub, b_ub=bub, bounds=(0, None), method="highs")
    hi = linprog(-np.ones_like(x), A_ub=Aub, b_ub=bub, bounds=(0, None), method="highs") if gamma is not None else None
    return (lo.fun if lo.status == 0 else np.nan,
            (-hi.fun if hi is not None and hi.status == 0 else np.inf))


def stieltjes_residual(N):
    """min max_n |int dmu/(lam+s_n) - E(n)| over mu >= 0 ; s_n = 1.5 n (physical units)."""
    lam = np.logspace(-4, 5, 3000)
    s = 1.5 * np.arange(1, N + 1)
    A = 1.0 / (lam[None, :] + s[:, None])
    K = len(lam)
    # variables mu (K), r (1): minimise r s.t. |A mu - E| <= r
    c = np.concatenate([np.zeros(K), [1.0]])
    Aub = np.vstack([np.hstack([A, -np.ones((N, 1))]), np.hstack([-A, -np.ones((N, 1))])])
    bub = np.concatenate([E[:N], -E[:N]])
    res = linprog(c, A_ub=Aub, b_ub=bub, bounds=(0, None), method="highs",
                  options=dict(primal_feasibility_tolerance=1e-10))
    return res.fun


if __name__ == "__main__":
    out = ["(H) Phi >= 0: rigorous lower bound on 4 pi eta/s from E(1..N) (true Einstein value 1)"]
    for N in range(1, 9):
        lo, _ = hausdorff_bounds(N)
        out.append(f"  N={N}: 4pi eta/s >= {lo:.4f}")
    out.append("(H) with Euclidean data known only to +-eps (lower bound, no M needed)")
    for eps in (1e-2, 1e-3, 1e-4, 1e-6):
        row = "  ".join(f"N={N}: {hausdorff_bounds(N, tol=eps)[0]:.4f}" for N in (2, 4, 6, 8))
        out.append(f"  eps={eps:.0e}: {row}")
    out.append("(H)+(M): two-sided bounds with int x^{-gamma} dnu <= M (gamma in units of y_1 = 1.5, i.e. Re s > -gamma)")
    for gamma, M in ((0.5, 1.5), (0.5, 2.0), (0.5, 5.0), (1.0, 5.0), (1.0, 20.0)):
        for N in (2, 4, 8):
            lo, hi = hausdorff_bounds(N, gamma, M)
            out.append(f"  gamma={gamma} M={M} N={N}: {lo:.4f} <= 4 pi eta/s <= {hi:.4f}")
    out.append("(S) Phi completely monotone (E Stieltjes): minimal sup-residual of a positive-mu fit to E(1..N)")
    for N in range(2, 9):
        out.append(f"  N={N}: residual = {stieltjes_residual(N):.3e}")
    txt = "\n".join(out)
    print(txt)
    open(__file__.replace("p1_einstein_positivity.py", "out_p1_einstein_positivity.txt"), "w").write(txt + "\n")
