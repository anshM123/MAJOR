"""Numerical falsification test of the Separation Lemma and of the acausality of pleated slices.

C = acausal circle in the boundary of AdS_3: points p(u) = (u,1)^T (phi(u),1), u in RP^1, phi an increasing homeomorphism
(here: a random increasing piecewise-Moebius-free map given by sorted random values).  Positive lifts: for
theta in (0,pi), a = (cos th, sin th), b = (cos ph, sin ph) with ph = Phi(th) increasing; then
N(p,p') = sin(th'-th) sin(ph'-ph) > 0 for th < th'.

Separation Lemma: cut C at two points into closed arcs alpha, beta.  For X a positive combination of points of alpha and
Y a positive combination of points of beta (both of positive determinant),
        tr(adj(X) Y) / (2 sqrt(det X det Y)) >= 1,
with equality only if X and Y are proportional and supported on the two cut points.
(>1 means the points x = X/sqrt(det X), y = Y/sqrt(det Y) of AdS_3 are joined by a spacelike geodesic.)

Tests: (1) random positive combinations over complementary arcs; (2) random pairs of ideal triangles with non-crossing
vertex sets, and pairs of points in the same triangle; (3) control: triangles with crossing vertex sets."""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
def adj(x): return np.array([[x[1, 1], -x[0, 1]], [-x[1, 0], x[0, 0]]])
def sep(X, Y): return np.trace(adj(X) @ Y)/(2*np.sqrt(np.linalg.det(X)*np.linalg.det(Y)))
def circle(n):
    th = np.sort(rng.uniform(0.02, np.pi-0.02, n)); ph = np.sort(rng.uniform(0.02, np.pi-0.02, n))
    return [np.outer([np.cos(t), np.sin(t)], [np.cos(p), np.sin(p)]) for t, p in zip(th, ph)]
def combo(P, idx):
    lam = rng.exponential(size=len(idx))*rng.choice([1e-3, 1, 1e3], size=len(idx))
    return sum(l*P[i] for l, i in zip(lam, idx))
worst = {1: np.inf, 2: np.inf, "same": np.inf}; bad = {1: 0, 2: 0, "same": 0}; ctrl_time = 0; ctrl_tot = 0
for trial in range(40000):
    n = int(rng.integers(6, 14)); P = circle(n)
    # (1) complementary arcs: cut at indices i<j; alpha = i..j, beta = j..n-1,0..i (cut points belong to both)
    i, j = sorted(rng.choice(n, 2, replace=False))
    A = list(range(i, j+1)); B = list(range(j, n)) + list(range(0, i+1))
    ia = rng.choice(A, size=min(len(A), int(rng.integers(2, 5))), replace=False); ib = rng.choice(B, size=min(len(B), int(rng.integers(2, 5))), replace=False)
    X, Y = combo(P, ia), combo(P, ib)
    if np.linalg.det(X) > 1e-300 and np.linalg.det(Y) > 1e-300:
        s = sep(X, Y); worst[1] = min(worst[1], s); bad[1] += s < 1 - 1e-9
    # (2) two ideal triangles with non-crossing vertex sets: choose 6 indices a<b<c<d<e<f; T=(a,b,f), T'=(c,d,e)  [all of T' in arc (b,f)]
    idx = np.sort(rng.choice(n, 6, replace=False)); T, Tp = [idx[0], idx[1], idx[5]], [idx[2], idx[3], idx[4]]
    s = sep(combo(P, T), combo(P, Tp)); worst[2] = min(worst[2], s); bad[2] += s < 1 - 1e-9
    s = sep(combo(P, T), combo(P, T)); worst["same"] = min(worst["same"], s); bad["same"] += s < 1 - 1e-9
    # adjacent triangles sharing an edge: T=(a,b,c), T'=(a,c,d)
    s = sep(combo(P, [idx[0], idx[1], idx[2]]), combo(P, [idx[0], idx[2], idx[3]])); worst[2] = min(worst[2], s); bad[2] += s < 1 - 1e-9
    # (3) control: crossing triangles T=(a,c,e), T'=(b,d,f): search a few random pairs for a timelike-related one
    hit = False
    for rep in range(20):
        if sep(combo(P, [idx[0], idx[2], idx[4]]), combo(P, [idx[1], idx[3], idx[5]])) < 1: hit = True; break
    ctrl_time += hit; ctrl_tot += 1
print(f"(1) complementary arcs:        violations {bad[1]}/40000, smallest ratio {worst[1]:.12f}")
print(f"(2) non-crossing triangles:    violations {bad[2]}/80000, smallest ratio {worst[2]:.12f}")
print(f"    same triangle:             violations {bad['same']}/40000, smallest ratio {worst['same']:.12f}")
print(f"(3) control, crossing triangles: a timelike-related pair of points found in {ctrl_time}/{ctrl_tot} configurations")
