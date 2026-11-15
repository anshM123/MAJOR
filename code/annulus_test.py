"""Falsification test of Theorem 1 (multi-metric RT) on the hyperbolic cylinder (rotating-BTZ topology).

Surface: H^2/<z -> e^ell z>.  One ideal boundary circle carries n marked points (coordinates 0 <= s_1<...<s_n < ell);
the other end of the cylinder belongs to the purifier.  The n elementary arcs a_k=(s_k,s_{k+1}) are assigned to
parties 0..N-1 or to the purifier.

Cut classes for a region A: a non-crossing perfect matching of the endpoints of A, plus the choice of the face of the
chord diagram that contains the hole of the cylinder; if that face lies on the A-side the core geodesic is added.
Regulated lengths (horocycles invariant under the deck group; additive constants per endpoint dropped):
   arc i--j homotopic to the boundary segment of coordinate length d :  2*log(2*sinh(d/2))
   core geodesic: ell.
Multi-metric functional: F(A) = min_c sum_j w_j * length_j(c)  with the SAME class c in all metrics.
Control (not of Crofton type, expected to violate): min_c max_j length_j(c)."""
import numpy as np, itertools, sys
from functools import lru_cache
from contraction import catalogue, REV_MMI

@lru_cache(None)
def matchings(pts):
    """non-crossing perfect matchings of a tuple of cyclically ordered labels"""
    if not pts: return ((),)
    out = []; a = pts[0]
    for k in range(1, len(pts), 2):
        for m1 in matchings(pts[1:k]):
            for m2 in matchings(pts[k+1:]):
                out.append(((a, pts[k]),) + m1 + m2)
    return tuple(out)

def region_classes(Aarcs, n):
    """yield (chords, hole_arc_index, core_needed) ; arcs indexed 0..n-1, arc k = (point k, point k+1 mod n)"""
    inA = [k in Aarcs for k in range(n)]
    ends = tuple(k for k in range(n) if inA[k-1] != inA[k])      # point k separates arcs k-1 and k
    for M in matchings(ends):
        seen = set()
        for h in range(n):
            # face signature of arc h: for each chord, which side h is on
            sig = tuple((i <= h < j) for (i, j) in M)
            if sig in seen: continue
            seen.add(sig); yield M, sig, inA[h]

LOG_INV_EPS = 0.0   # log(1/eps): each arc end contributes this much (horoball cut-off); 0 = finite parts only

def lengths(M, sig, core, ell, s):
    tot = (ell if core else 0.0) + 2*LOG_INV_EPS*len(M)
    for (i, j), inside in zip(M, sig):
        d = s[j] - s[i]
        if inside: d = ell - d            # hole on the [s_i,s_j] side -> arc is homotopic to the complementary segment
        tot += 2*np.log(2*np.sinh(d/2))
    return tot

def F(Aarcs, n, metrics, w, mode):
    best = np.inf
    for M, sig, core in region_classes(frozenset(Aarcs), n):
        ls = [lengths(M, sig, core, ell, s) for ell, s in metrics]
        if mode == "sum": v = sum(wi*li for wi, li in zip(w, ls))
        elif mode == "max": v = max(ls)
        best = min(best, v)
    return best

def run(name, ineq, k, trials, rng, mode):
    lhs, rhs, parties = ineq; N = len(parties); worst = np.inf; viol = 0
    for t in range(trials):
        n = int(rng.integers(N+1, N+5))
        # party assignment: every party gets >=1 arc, at least one purifier arc
        while True:
            lab = rng.integers(-1, N, size=n)
            if set(range(N)) <= set(lab.tolist()) and -1 in lab: break
        metrics = []
        for j in range(k):
            ell = float(np.exp(rng.uniform(np.log(0.3), np.log(12))))
            s = np.sort(rng.uniform(0, ell, size=n)); metrics.append((ell, s))
        w = rng.uniform(0.2, 1.0, size=k)
        def Sreg(A): return F([a for a in range(n) if lab[a] in A], n, metrics, w, mode)
        val = sum(c*Sreg(A) for c, A in lhs) - sum(c*Sreg(A) for c, A in rhs)
        scale = sum(c*abs(Sreg(A)) for c, A in lhs) + 1.0
        worst = min(worst, val); viol += (val < -1e-9*scale)
    print(f"{name:10s} k={k} mode={mode:4s}: trials={trials}  violations={viol}  min(LHS-RHS)={worst:+.3e}", flush=True)

if __name__ == "__main__":
    import annulus_test as _self
    if len(sys.argv) > 3: _self.LOG_INV_EPS = float(sys.argv[3]); globals()["LOG_INV_EPS"] = float(sys.argv[3])
    print("log(1/eps) per arc end =", LOG_INV_EPS)
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0); T = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    cat = catalogue()
    for name in ["SSA", "MMI", "cyclic5", "cyclic7"]:
        for k in (1, 2, 3):
            run(name, cat[name], k, T if name != "cyclic7" else T//4, rng, "sum")
    print("--- controls (functionals NOT covered by the theorem; violations expected for some) ---")
    for name in ["MMI", "cyclic5"]:
        run(name, cat[name], 2, T, rng, "max")
    run("revMMI", REV_MMI, 2, T, rng, "sum")
