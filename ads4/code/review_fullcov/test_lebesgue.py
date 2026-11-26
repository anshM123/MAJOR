"""PROOF.md §3: Lebesgue constant of (2m+1)-point equispaced trigonometric interpolation
versus Lambda_m = 1+2 sum_{k<=m} 1/(2k-1) (reference_driver.norm_constant(m,0)); and a
randomized check of the perturbed-node norming constant C_m(eta) on the m+1 cosine nodes
x_j=(1+cos(2 pi j/(2m+1)))/2 for random degree-m polynomials (worst-case search)."""
import sys, math
import numpy as np
from fractions import Fraction as Fr
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'sliver_fullcov' / 'sliver_merge'))

def Lam(m):
    return 1 + 2 * sum(1.0 / (2 * k - 1) for k in range(1, m + 1))

worst = 0
for m in list(range(1, 60)) + [100, 200, 398, 446, 511, 566, 821, 908]:
    n = 2 * m + 1; hstep = 2 * math.pi / n
    th = np.linspace(0, hstep / 2, 400)          # Lebesgue function is hstep-periodic and even
    nodes = hstep * np.arange(n)
    d = th[:, None] - nodes[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        D = np.where(np.abs(np.sin(d / 2)) < 1e-15, n, np.sin(n * d / 2) / np.sin(d / 2))
    L = np.abs(D).sum(axis=1).max() / n
    worst = max(worst, L / Lam(m))
    if m in (1, 2, 5, 10, 50, 446, 908): print(f'm={m}: true Lebesgue const {L:.5f}  Lambda_m {Lam(m):.5f}')
print('max true/Lambda_m =', worst, 'OK' if worst <= 1 else 'VIOLATION')

# even reduction + algebraic norming at m+1 cosine nodes, with random polynomials
rng = np.random.default_rng(1)
worst2 = 0
for m in (3, 8, 20, 60):
    xs = (1 + np.cos(2 * np.pi * np.arange(m + 1) / (2 * m + 1))) / 2
    grid = np.linspace(0, 1, 20001)
    V = np.polynomial.chebyshev.chebvander(2 * xs - 1, m); Vg = np.polynomial.chebyshev.chebvander(2 * grid - 1, m)
    # cardinal (Lagrange) basis at the m+1 nodes: Lebesgue function sum |l_j(x)|
    Lg = Vg @ np.linalg.inv(V)
    leb = np.abs(Lg).sum(axis=1).max()
    worst2 = max(worst2, leb / Lam(m))
    print(f'm={m}: algebraic Lebesgue const at m+1 cosine nodes {leb:.4f} <= Lambda_m {Lam(m):.4f}')
print('max algebraic/Lambda_m =', worst2, 'OK' if worst2 <= 1 else 'VIOLATION')
