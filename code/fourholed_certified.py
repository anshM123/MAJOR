"""Falsification test of Theorem 1 on the four-boundary wormhole topology (pi_1 = F_3), with structural guarantees.

"Certified" here means: guaranteed by construction plus checks in ordinary double precision with printed margins.
It is NOT interval arithmetic.  Known biases: base geometries are kept only if a non-boundary class enters a minimal
cut, and trials whose enumeration exceeds 400000 words are discarded.  The 'max' control computes the sum over
components of max_j length_j (not max_j of the sum) and is only a sanity control.

What is guaranteed/checked in each trial:
 (1) DISCRETE + FAITHFUL.  The group is a classical Schottky group: six pairwise disjoint closed intervals
     I_1,I_1',I_2,I_2',I_3,I_3' of the real line (in this order) and Moebius maps g_i with
     g_i(exterior of I_i) = interior of I_i'.  By the ping-pong lemma <g_1,g_2,g_3> is free of rank 3 and discrete,
     with fundamental domain the complement of the six half-discs.  We check g_i(l_i)=r_i', g_i(r_i)=l_i' and the
     disjointness margin.  There are NO relations to satisfy; the fourth boundary class is d=(g_3 g_2 g_1)^{-1}.
     The quotient has four funnels (boundary cycles of the side pairing: g_1, g_2, g_3 and g_3g_2g_1) and
     Euler characteristic -2, hence is a four-holed sphere.
 (2) ALL HOMOTOPY CLASSES + HOMOLOGY CONDITION.  Every closed geodesic is a cyclically reduced word.  If the word is
     x_1...x_n, the geodesic crosses the fundamental domain n times, the k-th time from side(x_k^{-1}) to side(x_{k+1}),
     so  length >= sum_k dist(side(x_k^{-1}), side(x_{k+1})).  A branch-and-bound over ALL cyclically reduced words with
     this lower bound <= Theta therefore visits every closed geodesic (simple or not) of weighted length <= Theta.
     The mod-2 homology class of a word is its vector of exponent-sum parities v in Z_2^3; the cut classes for a set T
     of boundaries are the collections of closed curves with total class v_T = sum_{i in T} [boundary_i]
     ([a]=(0,0,1),[b]=(0,1,0),[c]=(1,0,0) in (g1,g2,g3) coordinates, [d]=(1,1,1)); non-embedded collections never beat
     embedded ones (PROOFS.md Lemma 6.8).  Hence  F(T) = min over decompositions v_T = v_1+...+v_m of sum g(v_i),
     g(v) = min over words of class v of sum_j w_j length_j(word)   -- the SAME word in all metrics.
     Completeness: Theta is chosen >= max_v g(v) found, so no unvisited word can lower any g(v).
 (3) TIME DEPENDENCE.  For k=2 the two Schottky groups are sampled independently: X_L != X_R, so the spacetime has no
     time-reflection symmetry and each horizon carries angular momentum proportional to |l_L - l_R| (reported)."""
import numpy as np, sys, itertools
sys.setrecursionlimit(10000)

WMIN, WMAX, GMIN, GMAX = 0.05, 5.0, 0.05, 0.6
def mob(l, r): return np.array([[2.0, -(l+r)], [0.0, r-l]])
def sample_schottky(rng):
    # wide gaps between I_i and I_i' (long boundaries are avoided by moderate widths), moderate necks
    # widths spread over two decades: a large interval paired with a small one gives a long horizon, while two large
    # intervals separated only by small ones give a short neck (essential curve) -- the regime where non-boundary
    # classes enter the minimal cuts
    widths = np.exp(rng.uniform(np.log(WMIN), np.log(WMAX), size=6)); gaps = np.exp(rng.uniform(np.log(GMIN), np.log(GMAX), size=6))
    pos = 0.0; I = []
    for wd, gp in zip(widths, gaps): I.append((pos, pos+wd)); pos += wd+gp
    gens = []; S = np.array([[0.0, -1.0], [1.0, 0.0]])
    for i in range(3):
        (l, r), (lp, rp) = I[2*i], I[2*i+1]; t = rng.normal()*0.3
        M = np.array([[np.cosh(t), np.sinh(t)], [np.sinh(t), np.cosh(t)]])
        g = np.linalg.inv(mob(lp, rp)) @ M @ S @ mob(l, r); gens.append(g/np.sqrt(np.linalg.det(g)))
    return I, gens
def act(g, x): return (g[0, 0]*x+g[0, 1])/(g[1, 0]*x+g[1, 1])
def certify(I, gens):
    """returns (disjointness margin, max endpoint-mapping error)"""
    margin = min(I[k+1][0]-I[k][1] for k in range(5)); err = 0.0
    for i in range(3):
        (l, r), (lp, rp) = I[2*i], I[2*i+1]
        err = max(err, abs(act(gens[i], r)-lp), abs(act(gens[i], l)-rp))
        mid_out = l - 1.0                                # a point outside I_i must land strictly inside I_i'
        y = act(gens[i], mid_out); assert lp < y < rp
    return margin, err
def side_dist(I):
    c = [(a+b)/2 for a, b in I]; r = [(b-a)/2 for a, b in I]; D = np.zeros((6, 6))
    for p in range(6):
        for q in range(6):
            if p != q: D[p, q] = np.arccosh(((c[p]-c[q])**2 - r[p]**2 - r[q]**2)/(2*r[p]*r[q]))
    return D
# letters: 0,1,2 = g1,g2,g3 ; 3,4,5 = inverses.  exit side of letter x: g_i -> I_i' (index 2i+1), g_i^-1 -> I_i (2i)
def exit_side(x): return 2*x+1 if x < 3 else 2*(x-3)
def inv_letter(x): return (x+3) % 6
def length_from_trace(t): return 2*np.arccosh(max(abs(t)/2, 1.0))

def entropies_from_g(g):
    """F(v) = cheapest multiset of classes summing to v (shortest path in the group Z_2^3)"""
    vs = list(itertools.product((0, 1), repeat=3)); F = {v: (0.0 if v == (0, 0, 0) else np.inf) for v in vs}
    for _ in range(8):
        for v in vs:
            for u, (val, _) in g.items():
                vv = tuple((a+b) % 2 for a, b in zip(v, u))
                if F[vv] + val < F[v]: F[v] = F[vv] + val
    return F

def trial(rng, k, mode_weights="half"):
    data = [sample_schottky(rng) for _ in range(k)]
    margins = [certify(I, gens) for I, gens in data]
    mats = [[g for g in gens] + [np.linalg.inv(g) for g in gens] for _, gens in data]
    D = [side_dist(I) for I, _ in data]
    w = np.full(k, 1.0/k) if mode_weights == "half" else rng.uniform(0.3, 1.0, size=k)
    return data, margins, mats, D, w

if __name__ == "__main__":
    seed = int(sys.argv[1]); ntr = int(sys.argv[2]); rng = np.random.default_rng(seed)
    # --- calibration of the word/geodesic convention: the lower bound must hold for every word (it fails for the wrong order)
    data, margins, mats, D, w = trial(rng, 1)
    for order in (+1, -1):
        bad = 0; tot = 0
        for n in range(1, 6):
            for word in itertools.product(range(6), repeat=n):
                if any(word[i+1] == inv_letter(word[i]) for i in range(n-1)) or word[0] == inv_letter(word[-1]): continue
                P = np.eye(2)
                for x in (word if order > 0 else word[::-1]): P = P @ mats[0][x]
                lb = sum(D[0][exit_side(inv_letter(word[i])), exit_side(word[(i+1) % n])] for i in range(n))
                tot += 1; bad += length_from_trace(np.trace(P)) < lb - 1e-9
        print(f"calibration, product order {order:+d}: words (length<=5) violating the lower bound: {bad}/{tot}")

    # --- main test
    V = {"a": (0, 0, 1), "b": (0, 1, 0), "c": (1, 0, 0)}
    def vec(T): return tuple(sum(V[t][i] for t in T) % 2 for i in range(3))
    bd = {(0, 0, 1), (0, 1, 0), (1, 0, 0), (1, 1, 1)}
    def letter_mats(gens): return [g for g in gens] + [np.linalg.inv(g) for g in gens]
    def run(mats, costs, wcost, w, mode, theta):
        """exhaustive branch-and-bound; costs = side-distance matrices of the metrics whose marking is the Schottky one"""
        kk = len(mats); cost = np.zeros((6, 6))
        for x in range(6):
            for y in range(6):
                cost[x, y] = sum(wc*D[exit_side(inv_letter(x)), exit_side(y)] for wc, D in zip(wcost, costs)) if y != inv_letter(x) else np.inf
        mincost = cost[np.isfinite(cost)].min(); best = {}; visited = [0]
        def rec(word, prods, partial):
            first, last = word[0], word[-1]
            if last != inv_letter(first) and partial + cost[last, first] <= theta:
                visited[0] += 1
                ls = [length_from_trace(np.trace(prods[j])) for j in range(kk)]
                val = sum(w[j]*ls[j] for j in range(kk)) if mode == "sum" else max(ls)
                par = tuple(sum(1 for x in word if x % 3 == i) % 2 for i in range(3))
                if par != (0, 0, 0) and (par not in best or val < best[par][0]): best[par] = (val, tuple(word))
            if visited[0] > 400000: raise RuntimeError("too many words")
            for y in range(6):
                if y == inv_letter(last): continue
                p2 = partial + cost[last, y]
                if p2 + mincost > theta: continue
                word.append(y); rec(word, [prods[j] @ mats[j][y] for j in range(kk)], p2); word.pop()
        for x in range(6): rec([x], [mats[j][x] for j in range(kk)], 0.0)
        return best, visited[0]
    def upper(mats, w, mode, n_max=4):
        ub = {}; kk = len(mats)
        for n in range(1, n_max+1):
            for word in itertools.product(range(6), repeat=n):
                if any(word[i+1] == inv_letter(word[i]) for i in range(n-1)) or word[0] == inv_letter(word[-1]): continue
                P = [np.eye(2) for _ in range(kk)]
                for x in word: P = [P[j] @ mats[j][x] for j in range(kk)]
                ls = [length_from_trace(np.trace(P[j])) for j in range(kk)]
                val = sum(w[j]*ls[j] for j in range(kk)) if mode == "sum" else max(ls)
                par = tuple(sum(1 for x in word if x % 3 == i) % 2 for i in range(3))
                if par != (0, 0, 0) and (par not in ub or val < ub[par]): ub[par] = val
        return ub
    def full(mats, costs, wcost, w, mode):
        theta = max(upper(mats, w, mode).values())*(1+1e-9)
        best, nv = run(mats, costs, wcost, w, mode, theta)
        complete = (len(best) == 7) and max(v[0] for v in best.values()) <= theta
        F = entropies_from_g(best); Fb = entropies_from_g({p: v for p, v in best.items() if p in bd})
        return F, any(F[v] < Fb[v] - 1e-9 for v in F), nv, complete
    def jitter(I, rng, amp):
        pts = np.array([e for iv in I for e in iv]); gaps = np.diff(pts)
        gaps = gaps*np.exp(rng.normal(size=gaps.size)*amp); pts2 = pts[0] + np.concatenate([[0], np.cumsum(gaps)])
        return [(pts2[2*i], pts2[2*i+1]) for i in range(6)]
    def gens_from(I, rng):
        gens = []; S = np.array([[0.0, -1.0], [1.0, 0.0]])
        for i in range(3):
            (l, r), (lp, rp) = I[2*i], I[2*i+1]; t = rng.normal()*0.3
            M = np.array([[np.cosh(t), np.sinh(t)], [np.sinh(t), np.cosh(t)]])
            g = np.linalg.inv(mob(lp, rp)) @ M @ S @ mob(l, r); gens.append(g/np.sqrt(np.linalg.det(g)))
        return gens
    def twist(gens, which, sgn):
        """Dehn twist of the marking about g3g2 (which=0) or g2g1 (which=1); returns the images of g1,g2,g3"""
        g1, g2, g3 = gens; e = (g3 @ g2) if which == 0 else (g2 @ g1)
        if sgn < 0: e = np.linalg.inv(e)
        ei = np.linalg.inv(e); c = lambda g: e @ g @ ei
        return [g1, c(g2), c(g3)] if which == 0 else [c(g1), c(g2), g3]
    for k, mode in ((1, "sum"), (2, "sum"), (3, "sum"), (2, "max")):
        viol = 0; worst = np.inf; minmargin = np.inf; maxerr = 0; nwords = []; ess = 0; rot = []; incomplete = 0; done = 0; twisted = 0
        while done < ntr:
            I1, gens1 = sample_schottky(rng); m1 = letter_mats(gens1); D1 = side_dist(I1)
            try:
                F1, ess1, _, comp1 = full([m1], [D1], [1.0], [1.0], "sum")
            except RuntimeError: continue
            if not ess1: continue                      # keep only base geometries in which a neck matters
            mats = [m1]; marg = [certify(I1, gens1)]
            for j in range(1, k):
                Ij = jitter(I1, rng, 0.25); gj = gens_from(Ij, rng); marg.append(certify(Ij, gj))
                if rng.random() < 0.6:
                    gj = twist(gj, int(rng.integers(0, 2)), int(rng.choice([-1, 1]))); twisted += 1
                mats.append(letter_mats(gj))
            w = np.full(k, 1.0/k) if k == 2 else rng.uniform(0.3, 1.0, size=k)
            try:
                F, e_, nv, comp = full(mats, [D1], [w[0] if mode == "sum" else 1.0/k], w, mode)
            except RuntimeError: continue
            done += 1; nwords.append(nv); ess += e_; incomplete += (not comp)
            minmargin = min(minmargin, min(m for m, _ in marg)); maxerr = max(maxerr, max(e for _, e in marg))
            S = lambda *T: F[vec(T)]
            mmi = S("a", "b") + S("b", "c") + S("a", "c") - S("a") - S("b") - S("c") - S("a", "b", "c")
            ssa = min(S("a", "b") + S("b", "c") - S("b") - S("a", "b", "c"), S("a", "b") + S("a", "c") - S("a") - S("a", "b", "c"),
                      S("a", "c") + S("b", "c") - S("c") - S("a", "b", "c"))
            sc = sum(F.values()) + 1; worst = min(worst, mmi/sc, ssa/sc); viol += min(mmi, ssa) < -1e-9*sc
            if k == 2:
                for x in range(3): rot.append(abs(length_from_trace(np.trace(mats[0][x])) - length_from_trace(np.trace(mats[1][x]))))
        extra = f"; mean |l_L - l_R| of horizons = {np.mean(rot):.2f}" if rot else ""
        print(f"k={k} mode={mode}: trials={ntr}; violations of MMI/SSA = {viol}; min normalised (LHS-RHS) = {worst:+.1e}; "
              f"ping-pong margin >= {minmargin:.3f}, endpoint error <= {maxerr:.1e}; words visited: median {int(np.median(nwords))}, max {max(nwords)}; "
              f"incomplete enumerations = {incomplete}; trials where a non-boundary class enters a minimal cut = {ess}; "
              f"metrics with a relative Dehn twist = {twisted}{extra}", flush=True)
