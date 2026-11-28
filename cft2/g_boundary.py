"""Canonical finite part G versus boundary (CFT) expressions.
 G_can(A) = 4G S^eps(A) - sum_{p in dA} t_p^can(eps),  t_p^can from the universal-cover (lifted) neighbours of p
            [= (1/2)(l(p-,p)+l(p,p+)-l(p-,p+)) with lifted neighbours; equals log(1/(eps d_p)) of Prop. S4.3]
 G_CFT(A) = 4G S^eps(A) - sum_{p in dA} (1/2) 4G I(alpha_p^- : alpha_p^+),  alpha_p^pm = elementary intervals at p
Settings: (D) global AdS3 + boundary gravitons (vacuum Virasoro module), n>=3 points on the circle;
          (B) rotating BTZ + boundary-graviton-deformed slice, one boundary with n>=3 points + whole 2nd boundary;
          (F) free Dirac fermion vacuum on the circle (c=1).
Checks: G independent of eps and of the frame (random horocycle/Weyl normalisation c_i);
        G_CFT == G_can on (D); on (B) G_CFT != G_can when adjacent intervals are large; cone facets
        (all SA and all MMI instances, incl. purifier) for N=3,4 on G_can, G_CFT, S^eps.
Units: lengths l = 4 G_N S (chiral average).  c = 3/(2G_N) so 4G_N S = (6/c) S."""
import itertools, math, sys
import numpy as np

def matchings(pts):
    if not pts: yield []; return
    a = pts[0]
    for k in range(1, len(pts), 2):
        for m1 in matchings(pts[1:k]):
            for m2 in matchings(pts[k+1:]):
                yield [(a, pts[k])] + m1 + m2

def ivl(n, i, j):            # pieces covered going ccw from point i to point j
    out = []; k = i
    while k != j: out.append(k); k = (k+1) % n
    return frozenset(out)

class Geometry:
    """n marked points on one boundary circle. kind 'disc' (trivial holonomy) or 'btz' (hyperbolic, chiral
    horizon lengths a, b). arc(i,j) = regulated length of the geodesic homotopic to the ccw interval i->j."""
    def __init__(s, rng, n, kind):
        s.n, s.kind = n, kind
        s.cL = rng.uniform(0.2, 5, n); s.cR = rng.uniform(0.2, 5, n)     # frame / Weyl normalisations at each point
        if kind == 'disc':
            s.u = np.sort(rng.uniform(0, 2*math.pi, n)); s.v = np.sort(rng.uniform(0, 2*math.pi, n))
            s.v = (s.v - s.v[0] + s.u[0])            # same cyclic order (acausal slice), arbitrary positions
        else:
            s.a = rng.uniform(1, 7); s.b = rng.uniform(1, 7)
            s.th = np.sort(rng.uniform(0, s.a, n))
            amp = 0.8*s.b/(2*math.pi); ph = rng.uniform(0, 2*math.pi)
            s.phi = lambda t: (s.b/s.a)*t + amp*math.sin(2*math.pi*t/s.a + ph)   # graviton-deformed slice
    # chiral regulated length between lifted points (coordinate, normalisation)
    def chiral(s, side, ti, tj, ci, cj, eps):
        if s.kind == 'disc':
            return math.log((2*math.sin((tj-ti)/2))**2/(eps*ci*eps*cj))
        f = (lambda t: t) if side == 'L' else s.phi
        X1, X2 = math.exp(f(ti)), math.exp(f(tj))
        return math.log((X2-X1)**2/(eps*X1*ci*eps*X2*cj))
    def lift(s, side, i, shift):   # coordinate of point i lifted by 'shift' periods
        if s.kind == 'disc':
            t = s.u[i] if side == 'L' else s.v[i]; return t + 2*math.pi*shift
        return s.th[i] + s.a*shift
    def length(s, i, j, si, sj, eps):
        tot = 0
        for side, c in (('L', s.cL), ('R', s.cR)):
            tot += s.chiral(side, s.lift(side, i, si), s.lift(side, j, sj), c[i], c[j], eps)
        return tot/2
    def arc(s, i, j, eps):
        sj = 0 if (j > i) else 1
        return s.length(i, j, 0, sj, eps)
    def horizon(s): return 0.0 if s.kind == 'disc' else (s.a + s.b)/2
    def t_can(s, p, eps):
        n = s.n
        im, ip = (p-1) % n, (p+1) % n
        sm = -1 if p == 0 else 0; sp_ = 1 if p == n-1 else 0
        return 0.5*(s.length(im, p, sm, 0, eps) + s.length(p, ip, 0, sp_, eps) - s.length(im, ip, sm, sp_, eps))

def F(s, inA, b2, eps):
    """4G S(A): min over common embedded cut classes (non-crossing matchings, each arc homotopic to one side,
    nested-or-disjoint, parity colouring = A; horizon iff colouring near horizon differs from b2)."""
    n = s.n; E = [p for p in range(n) if inA[(p-1) % n] != inA[p]]
    best = math.inf
    for M in matchings(E):
        for orient in itertools.product((0, 1), repeat=len(M)):
            ivs = []; L = 0.0
            for (p, q), o in zip(M, orient):
                i, j = (p, q) if o == 0 else (q, p)
                ivs.append(ivl(n, i, j)); L += s.arc(i, j, eps)
            if any(not (u <= v or v <= u or not (u & v)) for u, v in itertools.combinations(ivs, 2)): continue
            for c0 in (0, 1):
                col = [(c0 + sum(1 for I in ivs if k in I)) % 2 for k in range(n)]
                if col != list(inA): continue
                h = (c0 + b2) % 2
                if s.kind == 'disc' and h: continue        # no second boundary / horizon in the disc
                best = min(best, L + h*s.horizon())
    return best

def entropies(s, lab, lab2, N, eps):
    regs = [frozenset(c) for r in range(1, N+1) for c in itertools.combinations(range(1, N+1), r)]
    n = s.n; out = {}
    tc = {p: s.t_can(p, eps) for p in range(n)}
    piece = lambda ks: [1 if k in ks else 0 for k in range(n)]
    tcft = {p: 0.5*(F(s, piece({(p-1) % n}), 0, eps) + F(s, piece({p}), 0, eps) - F(s, piece({(p-1) % n, p}), 0, eps))
            for p in range(n)}
    for A in regs:
        inA = [1 if lab[k] in A else 0 for k in range(n)]; b2 = 1 if lab2 in A else 0
        if s.kind == 'disc': b2 = 0
        f = F(s, inA, b2, eps); E = [p for p in range(n) if inA[(p-1) % n] != inA[p]]
        out[A] = (f, f - sum(tc[p] for p in E), f - sum(tcft[p] for p in E))
    return out, tc, tcft

def S_of(vals, X, N, k):     # entropy of a subset X of {0..N} (0 = purifier), pure total state
    X = frozenset(X)
    if 0 in X: X = frozenset(range(1, N+1)) - (X - {0})
    if not X: return 0.0
    return vals[X][k]

def instances(N):
    P = range(0, N+1); out = []
    subsets = [frozenset(c) for r in range(1, N+1) for c in itertools.combinations(P, r)]
    for X, Y in itertools.combinations(subsets, 2):
        if X & Y or len(X | Y) == N+1: continue
        out.append(('SA', [(X, 1), (Y, 1), (X | Y, -1)]))
    for X, Y, Z in itertools.combinations(subsets, 3):
        if X & Y or X & Z or Y & Z: continue
        out.append(('MMI', [(X | Y, 1), (Y | Z, 1), (X | Z, 1), (X, -1), (Y, -1), (Z, -1), (X | Y | Z, -1)]))
    return out

def labels(rng, n, N):      # cyclically adjacent pieces get different parties (P = entangling points only)
    while True:
        lab = [int(rng.integers(0, N+1)) for _ in range(n)]
        if all(lab[k] != lab[(k+1) % n] for k in range(n)): return lab

def run(kind, N, trials, seed):
    rng = np.random.default_rng(seed); inst = instances(N)
    worst = {0: math.inf, 1: math.inf, 2: math.inf}; drift = 0; frame = 0; gap = 0; nviol = 0; nbad_pt = 0
    for t in range(trials):
        n = int(rng.integers(3, 7)); s = Geometry(rng, n, kind)
        lab = labels(rng, n, N); lab2 = int(rng.integers(0, N+1))
        v1, tc1, tf1 = entropies(s, lab, lab2, N, 1e-6)
        v2, _, _ = entropies(s, lab, lab2, N, 1e-3)
        s.cL = s.cL*np.exp(rng.normal(0, 1, n)); s.cR = s.cR*np.exp(rng.normal(0, 1, n))   # change conformal frame
        v3, _, _ = entropies(s, lab, lab2, N, 1e-6)
        for A in v1:
            drift = max(drift, abs(v1[A][1]-v2[A][1]), abs(v1[A][2]-v2[A][2]))
            frame = max(frame, abs(v1[A][1]-v3[A][1]), abs(v1[A][2]-v3[A][2]))
            gap = max(gap, abs(v1[A][1]-v1[A][2]))
        nbad_pt += sum(1 for p in range(n) if abs(tc1[p]-tf1[p]) > 1e-9)
        scale = max(abs(v1[A][0]) for A in v1) + 1e-12
        bad = False
        for name, terms in inst:
            for k in (0, 1, 2):
                q = sum(c*S_of(v1, X, N, k) for X, c in terms)/scale
                worst[k] = min(worst[k], q)
                if k == 2 and q < -1e-9: bad = True
        nviol += bad
    print(f"[{kind}] N={N}, {trials} configs: eps-drift of G = {drift:.1e}; frame(Weyl)-change of G = {frame:.1e}; "
          f"max|G_can-G_CFT| = {gap:.3g}; points with t_CFT != t_can: {nbad_pt}")
    print(f"      min normalised SA/MMI value: S^eps {worst[0]:.2e}, G_can {worst[1]:.2e}, G_CFT {worst[2]:.2e}; "
          f"configs where G_CFT leaves the cone: {nviol}/{trials}")

def fermion(trials, seed):
    """Free Dirac fermion (c=1) vacuum on the circle: 6 S(A) = 2[sum log ch(a_i,b_j) - sum_{i<j} log ch(a_i,a_j)
       - sum_{i<j} log ch(b_i,b_j) - n log eps], ch = chord length.  (same units as above: 4G S <-> (6/c) S)"""
    rng = np.random.default_rng(seed); worst = math.inf; maxI3 = 0; worstpair = math.inf
    for t in range(trials):
        n = int(rng.integers(3, 8)); x = np.sort(rng.uniform(0, 2*math.pi, n)); N = 4
        lab = labels(rng, n, N); eps = 1e-6
        ch = lambda p, q: abs(2*math.sin((x[p]-x[q])/2))
        def S(inA):
            E = [p for p in range(n) if inA[(p-1) % n] != inA[p]]
            a = [p for p in E if inA[p] == 1]; b = [p for p in E if inA[p] == 0]   # left / right ends of A-intervals
            val = sum(math.log(ch(i, j)) for i in a for j in b) - sum(math.log(ch(i, j)) for i, j in itertools.combinations(a, 2)) \
                  - sum(math.log(ch(i, j)) for i, j in itertools.combinations(b, 2)) - len(a)*math.log(eps)
            return 2*val
        piece = lambda ks: [1 if k in ks else 0 for k in range(n)]
        I = lambda K1, K2: S(piece(K1)) + S(piece(K2)) - S(piece(K1 | K2))
        tcft = {p: 0.5*I({(p-1) % n}, {p}) for p in range(n)}
        vals = {}
        for r in range(1, N+1):
            for A in itertools.combinations(range(1, N+1), r):
                A = frozenset(A); inA = [1 if lab[k] in A else 0 for k in range(n)]
                E = [p for p in range(n) if inA[(p-1) % n] != inA[p]]
                G = S(inA) - sum(tcft[p] for p in E)
                # Bell-pair-sum prediction: sum over non-adjacent pieces alpha in A, beta not in A of (1/2) I(alpha:beta)
                pred = sum(0.5*I({al}, {be}) for al in range(n) for be in range(n)
                           if inA[al] and not inA[be] and (be-al) % n not in (1, n-1))
                worstpair = min(worstpair, -abs(G - pred))
                vals[A] = (G, G, G)
        for k1, k2, k3 in itertools.combinations(range(n), 3):
            maxI3 = max(maxI3, abs(I({k1}, {k2}) + I({k1}, {k3}) - I({k1}, {k2, k3})))
        scale = max(1.0, max(abs(v[0]) for v in vals.values()))   # n=3: G == 0 identically
        for name, terms in instances(N):
            worst = min(worst, sum(c*S_of(vals, X, N, 0) for X, c in terms)/scale)
    print(f"[fermion] {trials} configs, N=4: max|I3| over piece triples = {maxI3:.1e}; "
          f"max|G - Bell-pair-sum formula| = {-worstpair:.1e}; min normalised SA/MMI value of G = {worst:.2e}")

if __name__ == '__main__':
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    run('disc', 3, T, 1); run('disc', 4, T//2, 2)
    run('btz', 3, T, 3); run('btz', 4, T//2, 4)
    fermion(T, 5)
