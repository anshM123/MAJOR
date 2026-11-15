r"""Numerical test of Steps 3-4 of ads4/LIGHTCONE_THEOREM.md on synthetic positive measures.
Measure: towers at x = 2m (m = 0..5) made of many atoms with random spread and prescribed O(eps) shifts, an 'identity
image' atom at x = -2*Delta of weight eps*small, and junk atoms anywhere in [-1.5, 12].  We compute L(s) = \int e^{-s x} dM
exactly, apply the functionals K_1, K_2, and compare  \int_{|x|<eta} x e^{-S x} dM  (computed directly) with
-(1/(2P)) * [ \int K_1 dM ]  and with the planted first moment."""
import numpy as np
rng = np.random.default_rng(1)
def run(eps, M=4, kappa=0.3, eta=1.0, junk=1.0):
    S = 0.5*kappa*np.log(1/eps)                          # y = e^{-2S} = eps^kappa, as in the theorem (y = z^beta, eps = z^{tau_*/2}); y^M = eps^{1.2} = o(eps)
    xs, ws = [], []
    for m in range(6):                                   # towers
        n = 400; gam = -rng.uniform(0.5, 2.0, n)*eps      # anomalous dimensions, order eps, negative on average
        xs += list(2*m + gam); ws += list(rng.uniform(0.5, 1.5, n)/n*np.exp(-0.3*m))
    xs += [-1.3]; ws += [junk*eps**1.8]                   # identity image, o(eps)
    k = 200; xs += list(rng.uniform(-1.5, 12, k)); ws += list(junk*eps**1.8*rng.uniform(0, 1, k)/k)   # junk, o(eps)
    xs, ws = np.array(xs), np.array(ws)
    Q = lambda x: np.prod([4*np.sinh((x-2*m)/2)**2 for m in range(1, M)], axis=0)
    P = Q(0.0)
    K1 = np.exp(-S*xs)*2*np.sinh(xs)*Q(xs); K2 = np.exp(-S*xs)*4*np.sinh(xs/2)**2*Q(xs)
    win = np.abs(xs) < eta
    direct = np.sum(ws[win]*xs[win]*np.exp(-S*xs[win])); via = np.sum(ws*K1)/(2*P)
    return direct, via, np.sum(ws*K2), np.sum(ws[win]*xs[win]**2*np.exp(-S*xs[win]))
for label, junk in (("with o(eps) mass at negative x and between towers (weights eps^1.8)", 1.0), ("towers only (no junk, no identity image)", 0.0)):
    print(label)
    print("   eps       direct first moment     via K_1 functional     relative difference     int K_2 dM")
    for eps in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14):
        d, v, k2, m2 = run(eps, junk=junk)
        print(f"{eps:8.0e}   {d: .8e}       {v: .8e}       {abs(d-v)/abs(d):.3e}          {k2:.3e}")
print("Expected: relative difference -> 0 as eps -> 0 (rate eps^0.2 from the towers m >= M, weight y^M = eps^1.2), with large")
print("constants, because the kernels K_1, K_2 are large at negative x and at the uncancelled towers.")

print("\nM = 6 (every planted tower is a double zero of the kernels), towers only:")
for eps in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
    d, v, k2, m2 = run(eps, M=6, junk=0.0); print(f"  eps={eps:.0e}  direct={d:.6e}  via K_1={v:.6e}  rel.diff={abs(d-v)/abs(d):.2e}")
print("M = 6, with o(eps) junk (weights eps^1.8) at negative x and between towers:")
for eps in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10):
    d, v, k2, m2 = run(eps, M=6, junk=1.0); print(f"  eps={eps:.0e}  direct={d:.6e}  via K_1={v:.6e}  rel.diff={abs(d-v)/abs(d):.2e}")
