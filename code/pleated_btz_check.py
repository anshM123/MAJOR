"""Numerical test of Proposition P (S3) on a rotating BTZ geometry, including triangles spinning around the horizon.

Group: h = (u -> e^a u, v -> e^b v), a != b (rotating).  Acausal invariant circle C: v = phi(u), phi(u) = sign(u)|u|^{b/a},
with p(u) = (u,1)^T (phi(u),1), p(inf) = E11.  The horizon geodesic g(A) joins p(0) and p(inf).
Pleated slice pieces tested:
  side 1 (u>0): fan triangles (e1, v_k, v_{k+1}),  v_k = e^{a k} v_0,  e1 in {0, inf}  (spinning direction)
                 + dyadic triangles in the half-plane over [v_k, v_{k+1}]
  side 2 (u<0): fan triangles (e2, w_k, w_{k+1}),  w_k = -e^{a k} w_0, e2 in {0, inf}
                 + dyadic triangles
  the axis itself.
For random pairs of points on distinct pieces (and on the same piece) we evaluate
  ratio = tr(adj(X) Y) / (2 sqrt(det X det Y))   (>1 <=> spacelike separated; the lemma predicts >=1, =1 only at equal points)."""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
def P(u, a, b):
    if u == np.inf: return np.array([[1.0, 0], [0, 0]])
    ph = np.sign(u)*abs(u)**(b/a); return np.array([[u*ph, u], [ph, 1.0]])
adj = lambda x: np.array([[x[1, 1], -x[0, 1]], [-x[1, 0], x[0, 0]]])
def ratio(X, Y): return np.trace(adj(X) @ Y)/(2*np.sqrt(np.linalg.det(X)*np.linalg.det(Y)))
def point(verts, a, b):
    lam = rng.exponential(size=len(verts))*rng.choice([1e-2, 1, 1e2], size=len(verts))
    X = sum(l*P(u, a, b)/np.abs(P(u, a, b)).max() for l, u in zip(lam, verts)); return X
def dyadic(lo, hi, depth):
    out = []; 
    def rec(l, r, d):
        if d == 0: return
        m = np.sign(l)*np.sqrt(l*r) if l*r > 0 else (l+r)/2
        out.append((l, m, r)); rec(l, m, d-1); rec(m, r, d-1)
    rec(lo, hi, depth); return out
worst = np.inf; bad = 0; tot = 0
for trial in range(300):
    a = rng.uniform(0.4, 2.0); b = rng.uniform(0.4, 2.0); v0 = np.exp(rng.normal()); w0 = np.exp(rng.normal())
    e1 = rng.choice([0.0, np.inf]); e2 = rng.choice([0.0, np.inf]); K = 12
    pieces = [(0.0, np.inf)]                                           # the axis
    for k in range(-K, K):
        vk, vk1 = np.exp(a*k)*v0, np.exp(a*(k+1))*v0; wk, wk1 = -np.exp(a*k)*w0, -np.exp(a*(k+1))*w0
        pieces.append((e1, vk, vk1)); pieces.append((e2, wk1, wk))
        if abs(k) <= 2:
            pieces += dyadic(vk, vk1, 3) + [tuple(sorted(t)) for t in dyadic(wk1, wk, 3)]
    for rep in range(400):
        i, j = rng.integers(len(pieces), size=2)
        r = ratio(point(pieces[i], a, b), point(pieces[j], a, b)); tot += 1
        if not np.isfinite(r): continue
        worst = min(worst, r); bad += r < 1 - 1e-9
# control: a triangle that crosses the axis (vertices on both sides) against the fan triangles
cbad = ctot = 0
for trial in range(300):
    a = rng.uniform(0.4, 2.0); b = rng.uniform(0.4, 2.0)
    for rep in range(100):
        v = np.exp(rng.normal(size=4))
        r = ratio(point((-v[0], v[1], -v[2]*0 + v[2]), a, b), point((0.0, -v[3], -v[3]*np.exp(a)), a, b)); ctot += 1
        cbad += r < 1 - 1e-9
print(f"control (crossing triangles): {cbad}/{ctot} causally related pairs")
print(f"pairs tested: {tot}; pairs with ratio < 1 (causally related): {bad}; smallest ratio: {worst:.12f}")
