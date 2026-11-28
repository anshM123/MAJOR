"""Exact tests in the simply connected case (X = H^2, k metrics = k configurations of the same cyclically
ordered marked points; this is pure AdS3 with k=2, w=1/2 when the two chiral configurations are the
left/right boundary coordinates).  Tests, to machine precision:
 (R)/(D): F^eps(A) - sum_{p in dA} t_p(eps) == G(A) := min_c sum_tau nu(tau) i(tau,c), with nu(tau) the
          Liouville masses of the pairs of non-adjacent pieces (the classes tau_p removed);
 (3*):    G(A) == min-cut of the explicit cell graph of the arrangement of minimisers (Sec. 3 of FULLCONE.md),
          with nu-hat = sum_tau nu(tau) * (Dirac mass on one generic geodesic of class tau);
 cone:    F^eps (small eps) and G satisfy all facets of C_N, N=3,4,5.
Coordinates: upper half-plane boundary R; marked points x_0<...<x_{n-1}; piece k=(x_k,x_{k+1}), piece n-1
contains infinity.  Horoball at x_k of Euclidean diameter eps*s_k (random s_k: arbitrary normalisation)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, itertools, functools, fc_common as fc

@functools.lru_cache(None)
def matchings(pts):
    if not pts: return ((),)
    out = []; a = pts[0]
    for k in range(1, len(pts), 2):
        for m1 in matchings(pts[1:k]):
            for m2 in matchings(pts[k+1:]):
                out.append(((a, pts[k]),) + m1 + m2)
    return tuple(out)

def arc_len(x, h, i, j): return np.log((x[i]-x[j])**2/(h[i]*h[j]))

def liouville_pair(x, I, J):
    """m(I x J) for pieces I<J (indices), boundary points on R; piece n-1 contains infinity."""
    n = len(x)
    def ends(k): return (x[k], x[k+1]) if k < n-1 else (x[n-1], x[0])   # last piece wraps through infinity
    a1, a2 = ends(I); b1, b2 = ends(J)
    def cr(p, q, r, s):   # log of cross ratio for intervals (p,q),(r,s) in cyclic order p<q<r<s (with inf allowed via wrap)
        return np.log(abs((r-p)*(s-q)/((r-q)*(s-p))))
    return cr(a1, a2, b1, b2)

def d_p(x, k):
    """chart z -> -1/(z - x_k): q_pm = extremes of images of the other marked points"""
    q = [-1.0/(x[r]-x[k]) for r in range(len(x)) if r != k]
    return max(q) - min(q)

def setup(rng, N, n, kmet):
    while True:
        lab = rng.integers(-1, N, size=n)
        if set(range(N)) <= set(lab.tolist()) and -1 in lab: break
    lab = [N if l < 0 else int(l) for l in lab]
    mets = []
    for j in range(kmet):
        x = np.sort(rng.normal(0, 1, size=n)*np.exp(rng.uniform(-1, 1)))
        s = np.exp(rng.uniform(-1, 1, size=n))
        mets.append((x, s))
    w = rng.uniform(0.2, 1.0, size=kmet)
    return lab, mets, w

def endpoints(lab, n, m, N):
    inA = [(lab[k] < N and (m >> lab[k]) & 1) for k in range(n)]
    return tuple(k for k in range(n) if inA[k-1] != inA[k]), inA   # point k separates piece k-1 and piece k

def sep(i, j, I, J):
    """chord between points i<j separates pieces I,J? pieces between points i..j-1 are on one side"""
    inside = lambda K: i <= K < j
    return inside(I) != inside(J)

def trial(rng, N, n, kmet, eps=1e-6):
    lab, mets, w = setup(rng, N, n, kmet)
    pairs = [(I, J) for I in range(n) for J in range(I+1, n) if (J-I) % n not in (1, n-1)]
    nu = {pr: sum(wj*liouville_pair(x, *pr) for wj, (x, s) in zip(w, mets)) for pr in pairs}
    t = [sum(wj*np.log(1/(eps*s[k]*d_p(x, k))) for wj, (x, s) in zip(w, mets)) for k in range(n)]
    F, G, Gc, argmin = {}, {}, {}, {}
    for m in range(1, 2**N):
        ends, inA = endpoints(lab, n, m, N)
        best, bestc, bm = np.inf, np.inf, None
        for M in matchings(ends):
            v = sum(wj*sum(arc_len(x, eps*s, i, j) for i, j in M) for wj, (x, s) in zip(w, mets))
            vc = sum(nu[pr]*sum(sep(i, j, *pr) for i, j in M) for pr in pairs)
            if v < best: best, bm = v, M
            bestc = min(bestc, vc)
        F[m] = best; G[m] = best - sum(t[k] for k in ends); Gc[m] = bestc; argmin[m] = bm
    # cell graph of the arrangement of minimisers, in the Klein disc of metric 1
    x1 = mets[0][0]
    ang = lambda u: 2*np.arctan(u)
    P = np.array([[np.cos(ang(u)), np.sin(ang(u))] for u in x1])
    chords = sorted({c for M in argmin.values() for c in M})
    def side(c, z):
        a, b = P[c[0]], P[c[1]]
        return np.sign((b[0]-a[0])*(z[1]-a[1]) - (b[1]-a[1])*(z[0]-a[0]))
    def sv(z): return tuple(int(side(c, z)) for c in chords)
    def seg_inter(a, b, c, d):
        r = b-a; s = d-c; den = r[0]*s[1]-r[1]*s[0]
        if abs(den) < 1e-15: return None
        u = ((c-a)[0]*s[1]-(c-a)[1]*s[0])/den; v = ((c-a)[0]*r[1]-(c-a)[1]*r[0])/den
        return u if (0 < u < 1 and 0 < v < 1) else None
    # representative geodesic of each class: chord between generic points of the two pieces (metric 1)
    def piece_pt(K, frac):
        if K < n-1: u = x1[K] + frac*(x1[K+1]-x1[K])
        else:
            a0, a1 = ang(x1[n-1]), ang(x1[0]) + 2*np.pi
            th = a0 + frac*(a1-a0); return np.array([np.cos(th), np.sin(th)])
        return np.array([np.cos(ang(u)), np.sin(ang(u))])
    reps = {pr: (piece_pt(pr[0], rng.uniform(.3, .7)), piece_pt(pr[1], rng.uniform(.3, .7))) for pr in pairs}
    nodes, edges = set(), []
    for ci, c in enumerate(chords):
        a, b = P[c[0]], P[c[1]]
        us = sorted([0.0, 1.0] + [u for cj, d in enumerate(chords) if cj != ci
                                   for u in [seg_inter(a, b, P[d[0]], P[d[1]])] if u is not None])
        for u0, u1 in zip(us[:-1], us[1:]):
            mid = a + (u0+u1)/2*(b-a); nrm = np.array([-(b-a)[1], (b-a)[0]])*1e-7
            c1, c2 = sv(mid+nrm), sv(mid-nrm)
            e0, e1 = a+u0*(b-a), a+u1*(b-a)
            wt = sum(nu[pr]*(seg_inter(e0, e1, *reps[pr]) is not None) for pr in pairs)
            nodes |= {c1, c2}; edges.append((c1, c2, wt))
    label = {}
    for K in range(n):
        v = tuple(-1 if c[0] <= K < c[1] else 1 for c in chords); nodes.add(v)   # combinatorial side (ccw arc is on the right)
        if v in label and label[v] != lab[K]: raise RuntimeError(f"cell with two parties: K={K} v={v} lab={lab} chords={chords}")
        label[v] = lab[K]
    nodes = sorted(nodes)
    Gg = fc.graph_entropy(nodes, edges, label, N)
    eps1 = min(np.prod([(s[k]*d_p(x, k))**(-wj/w.sum()) for wj, (x, s) in zip(w, mets)]) for k in range(n))
    return dict(dRD=max(abs(G[m]-Gc[m]) for m in G), dGraph=max(abs(G[m]-Gg[m]) for m in G),
                tmin=min(t), F=F, G=G, eps1=eps1, lab=lab, mets=mets, w=w, n=n)

if __name__ == "__main__":
    rng = np.random.default_rng(7)
    for N, trials in ((3, 300), (4, 200), (5, 60)):
        stats = dict(dRD=0, dGraph=0, fF=np.inf, fG=np.inf, fFbig=np.inf)
        for _ in range(trials):
            n = int(rng.integers(N+2, N+5)); kmet = int(rng.integers(1, 4))
            r = trial(rng, N, n, kmet)
            stats["dRD"] = max(stats["dRD"], r["dRD"]); stats["dGraph"] = max(stats["dGraph"], r["dGraph"])
            stats["fF"] = min(stats["fF"], fc.facet_check(r["F"], N)[0]); stats["fG"] = min(stats["fG"], fc.facet_check(r["G"], N)[0])
        print(f"N={N} trials={trials}: max|G-G_classes|={stats['dRD']:.2e}  max|G-mincut(cell graph)|={stats['dGraph']:.2e}"
              f"  min facet(F^eps)={stats['fF']:+.2e}  min facet(G)={stats['fG']:+.2e}", flush=True)
