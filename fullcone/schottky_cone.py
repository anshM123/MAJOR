"""Full-cone test on multiply connected surfaces with interval parties (and whole-boundary parties).

Surfaces: rank-2 Schottky groups in the upper half-plane.  Four disjoint discs on R; generator g_k maps the
outside of D_k onto the inside of D_k' (g(z) = c' - r r'/(z - c), det = r r' > 0).  Disc order
D1 D1' D2 D2' gives a pair of pants (3 funnels), order D1 D2 D1' D2' a one-holed torus (1 funnel).
Several metrics share the combinatorics (same disc order, marked points on the same free arcs in the same
order), which fixes the marking; topological predicates are evaluated in metric 1 (and cross-checked in
metric 2).

For every region A the script computes, by exhaustive search over embedded cuts built from
  * arcs  (p_i, g p_j), g a reduced word of length <= L_ARC,
  * simple closed geodesics, primitive cyclically reduced words of length <= L_CL,
F^eps(A) = min_c sum_j w_j l^eps_{X_j}(c)   (same class in all metrics),
then t_p(eps) = sum_j w_j log(1/(h_p^(j) d_p^(j))) (FULLCONE.md Prop. 1) and G = F^eps - sum_{p in dA} t_p,
and checks: G independent of eps; all facets of C_N for F^eps (small eps) and for G; stability under
L_ARC -> L_ARC+1.  Controls: the 'max' functional min_c max_j l_j, and G shifted by -0.5*W*sum_p b_p."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, itertools, fc_common as fc

INV = [1, 0, 3, 2]   # letters: 0=a 1=A 2=b 3=B

def words(L):
    out = [()]
    frontier = [()]
    for _ in range(L):
        new = []
        for w in frontier:
            for k in range(4):
                if w and INV[w[-1]] == k: continue
                new.append(w + (k,))
        out += new; frontier = new
    return out

def inv_word(w): return tuple(INV[k] for k in reversed(w))

def gen_matrix(c, r, c2, r2):
    M = np.array([[c2, -r*r2 - c*c2], [1.0, -c]]); return M/np.sqrt(np.linalg.det(M))

class Surface:
    def __init__(self, discs, order):
        """discs: dict name->(c,r) for 'a','A','b','B' meaning D1, D1', D2, D2'. g_a: outside D_a -> inside D_A."""
        self.discs = discs
        ga = gen_matrix(*discs['a'], *discs['A']); gb = gen_matrix(*discs['b'], *discs['B'])
        self.gens = [ga, np.linalg.inv(ga), gb, np.linalg.inv(gb)]
        self.order = order
    def mat(self, w):
        M = np.eye(2)
        for k in w: M = M @ self.gens[k]
        return M
    @staticmethod
    def act(M, x): return (M[0, 0]*x + M[0, 1])/(M[1, 0]*x + M[1, 1])
    @staticmethod
    def deriv(M, x): return 1.0/(M[1, 0]*x + M[1, 1])**2
    def free_arcs(self):
        """free arcs of R between consecutive discs (sorted); the last one passes through infinity"""
        iv = sorted((c - r, c + r, nm) for nm, (c, r) in self.discs.items())
        arcs = [(iv[k][1], iv[k+1][0], iv[k][2], iv[k+1][2]) for k in range(3)]
        arcs.append((iv[3][1], iv[0][0], iv[3][2], iv[0][2]))   # through infinity
        return arcs

def random_surface(rng, topo, jitter=None):
    """topo 'pants': order a A b B ; 'torus': a b A B. returns discs dict."""
    names = ['a', 'A', 'b', 'B'] if topo == 'pants' else ['a', 'b', 'A', 'B']
    while True:
        gaps = rng.uniform(0.1, 0.6, size=4); rad = rng.uniform(0.3, 0.9, size=4)
        x = 0.0; discs = {}
        for nm, gp, rr in zip(names, gaps, rad):
            x += gp + rr; discs[nm] = (x, rr); x += rr
        # centre the configuration
        sh = x/2; discs = {k: (c - sh, r) for k, (c, r) in discs.items()}
        return discs

def boundary_circles(S):
    """free arcs glued into ideal boundary circles. A free arc runs left->right from the right end of disc X
    to the left end of disc Y. Left end of D_k is mapped by g_k to the right end of D_k' (and conversely)."""
    arcs = S.free_arcs()
    # map: arc k ends at left point of disc Y; under the pairing the left point of Y is identified with the
    # right point of partner(Y), which is the start of the arc beginning at partner(Y)'s right end.
    partner = {'a': 'A', 'A': 'a', 'b': 'B', 'B': 'b'}
    start_at = {a[2]: k for k, a in enumerate(arcs)}       # arc starting at right end of disc
    nxt = {k: start_at[partner[a[3]]] for k, a in enumerate(arcs)}
    seen, circles = set(), []
    for k in range(4):
        if k in seen: continue
        cyc = []; j = k
        while j not in seen: seen.add(j); cyc.append(j); j = nxt[j]
        circles.append(cyc)
    return arcs, circles

def arc_point(arcs, k, frac):
    lo, hi = arcs[k][0], arcs[k][1]
    if k == 3:   # through infinity: use the finite part to the right of lo (frac in (0,1) -> (lo, lo+span))
        return lo + frac*1.5
    return lo + frac*(hi - lo)

def linked(a, b, c, d, tol=1e-9):
    lo, hi = min(a, b), max(a, b)
    def inside(z): return lo + tol < z < hi - tol
    def onend(z): return abs(z - lo) < tol or abs(z - hi) < tol
    if onend(c) or onend(d): return False
    return inside(c) != inside(d)

class Config:
    def __init__(self, rng, topo, kmet, N, nmark, L_ARC=3, L_CL=3, L_CHK=8):
        base = random_surface(rng, topo)
        self.mets = []
        for j in range(kmet):
            if j == 0: discs = base
            else:   # independent geometry with the same combinatorics
                discs = random_surface(rng, topo)
            self.mets.append(Surface(discs, topo))
        S1 = self.mets[0]
        self.arcs, self.circles = boundary_circles(S1)
        self.topo, self.N = topo, N
        # marked points: choose circles that carry marked points; positions as fractions on free arcs
        ncirc = len(self.circles)
        while True:
            alloc = rng.integers(0, ncirc, size=nmark)
            marked_circles = sorted(set(alloc.tolist()))
            if topo == 'torus' or len(marked_circles) >= 1: break
        self.marks = []   # list of (free arc index, frac)
        for ci in range(ncirc):
            m = int((alloc == ci).sum())
            if m == 0: continue
            arcs_c = self.circles[ci]
            ks = rng.integers(0, len(arcs_c), size=m)
            for kk in ks: self.marks.append((arcs_c[kk], float(rng.uniform(0.1, 0.9))))
        self.marks.sort()
        # positions in each metric (same free arc, same fraction -> same cyclic order)
        self.pos = [[arc_point(boundary_circles(Sj)[0], k, f) for k, f in self.marks] for Sj in self.mets]
        self.hscale = [rng.uniform(0.5, 2.0, size=len(self.marks)) for _ in self.mets]
        # pieces: walk each circle, collect segments between consecutive marked points
        self.pieces = []; self.adj_pieces = {}
        for ci, cyc in enumerate(self.circles):
            L = len(cyc)
            pts = sorted((cyc.index(k), f, mi) for mi, (k, f) in enumerate(self.marks) if k in cyc)
            if not pts:
                self.pieces.append([(k, 0.0, 1.0) for k in cyc]); continue
            m = len(pts)
            for t in range(m):
                (t0, f0, i0), (t1, f1, i1) = pts[t], pts[(t+1) % m]
                if t0 == t1 and f1 > f0 and m > 1:
                    segs = [(cyc[t0], f0, f1)]
                else:
                    segs = [(cyc[t0], f0, 1.0)]; u = (t0 + 1) % L
                    while u != t1: segs.append((cyc[u], 0.0, 1.0)); u = (u + 1) % L
                    segs.append((cyc[t1], 0.0, f1))
                pid = len(self.pieces); self.pieces.append(segs)
                self.adj_pieces.setdefault(i0, [None, None])[1] = pid
                self.adj_pieces.setdefault(i1, [None, None])[0] = pid
        self.labels = None
        # words
        self.W_ARC = words(L_ARC); self.W_CHK = words(L_CHK); self.L_CL = L_CL
        self.mats = [[Sj.mat(w) for w in self.W_CHK] for Sj in self.mets]
        self.widx = {w: i for i, w in enumerate(self.W_CHK)}

    def assign(self, rng):
        """random party labels on pieces, every party and the purifier used"""
        P = len(self.pieces)
        if P < self.N + 1: return False
        while True:
            lab = rng.integers(0, self.N + 1, size=P)
            if set(range(self.N + 1)) <= set(lab.tolist()): break
        self.labels = [int(l) for l in lab]; return True

    # ---- geometry of components (metric j) ----
    def endpoints(self, comp, j):
        kind = comp[0]
        if kind == 'arc':
            _, i, w, k = comp
            M = self.mats[j][self.widx[w]]
            return self.pos[j][i], Surface.act(M, self.pos[j][k])
        _, w = comp
        M = self.mats[j][self.widx[w]]
        a, b, c, d = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
        disc = np.sqrt((a + d)**2 - 4)
        return ((a - d) + disc)/(2*c), ((a - d) - disc)/(2*c)

    def length(self, comp, j, eps):
        if comp[0] == 'arc':
            _, i, w, k = comp
            M = self.mats[j][self.widx[w]]
            u = self.pos[j][i]; v = Surface.act(M, self.pos[j][k])
            hu = eps*self.hscale[j][i]; hv = eps*self.hscale[j][k]*Surface.deriv(M, self.pos[j][k])
            return np.log((u - v)**2/(hu*hv))
        M = self.mats[j][self.widx[comp[1]]]
        return 2*np.arccosh(abs(M[0, 0] + M[1, 1])/2)

    def translates(self, comp, j, maxlen):
        u, v = self.endpoints(comp, j)
        out = []
        for w, M in zip(self.W_CHK, self.mats[j]):
            if len(w) > maxlen: continue
            out.append((Surface.act(M, u), Surface.act(M, v), w))
        return out

def build_components(C):
    j = 0
    comps = []
    n = len(C.marks)
    seen = set()
    for i in range(n):
        for k in range(n):
            for w in C.W_ARC:
                if not w and i == k: continue
                key = min((i, w, k), (k, inv_word(w), i))
                if key in seen: continue
                seen.add(key); comps.append(('arc',) + key)
    # closed: primitive cyclically reduced words up to L_CL, modulo rotation and inversion
    cseen = set()
    for w in words(C.L_CL):
        if not w or INV[w[0]] == w[-1]: continue
        L = len(w)
        if any(L % d == 0 and w == w[:d]*(L//d) for d in range(1, L)): continue
        rots = [w[t:] + w[:t] for t in range(L)]; iw = inv_word(w)
        rots += [iw[t:] + iw[:t] for t in range(L)]
        key = min(rots)
        if key in cseen: continue
        cseen.add(key); comps.append(('closed', key))
    return comps

def comp_len(comp): return len(comp[2]) if comp[0] == 'arc' else len(comp[1])

def linked_vec(u, v, A, B, tol=1e-11):
    lo, hi = min(u, v), max(u, v)
    sc = tol*max(1.0, abs(lo), abs(hi))
    onA = (np.abs(A - lo) < sc) | (np.abs(A - hi) < sc); onB = (np.abs(B - lo) < sc) | (np.abs(B - hi) < sc)
    inA = (A > lo + sc) & (A < hi - sc); inB = (B > lo + sc) & (B < hi - sc)
    return (inA != inB) & ~onA & ~onB

class Topology:
    """simplicity / disjointness / homology from lifted endpoints in metric j (exact up to the word bound:
    a translate h.c2 can link the base lift of c1 only if |h| <= |c1|+|c2|+2, see FULLCONE.md Sec. 7)"""
    def __init__(self, C, comps, j=0):
        self.C, self.j, self.comps = C, j, comps
        self.lens = np.array([len(w) for w in C.W_CHK])
        self.M = np.array(C.mats[j])
        self.tr = {}; self.cache = {}
    def base(self, ci): return self.C.endpoints(self.comps[ci], self.j)
    def lifts(self, ci):
        if ci not in self.tr:
            u, v = self.base(ci); M = self.M
            A = (M[:, 0, 0]*u + M[:, 0, 1])/(M[:, 1, 0]*u + M[:, 1, 1])
            B = (M[:, 0, 0]*v + M[:, 0, 1])/(M[:, 1, 0]*v + M[:, 1, 1])
            self.tr[ci] = (A, B)
        return self.tr[ci]
    def crosses(self, ci, cj):
        key = (min(ci, cj), max(ci, cj))
        if key in self.cache: return self.cache[key]
        bound = comp_len(self.comps[ci]) + comp_len(self.comps[cj]) + 2
        if bound > self.lens.max(): raise ValueError("word bound too small")
        sel = self.lens <= bound
        u, v = self.base(ci); A, B = self.lifts(cj)
        L = linked_vec(u, v, A[sel], B[sel])
        if ci == cj:   # remove translates that coincide with the base (stabiliser of a closed axis)
            sc = 1e-11*max(1.0, abs(u), abs(v))
            same = ((np.abs(A[sel]-u) < sc) & (np.abs(B[sel]-v) < sc)) | ((np.abs(A[sel]-v) < sc) & (np.abs(B[sel]-u) < sc))
            L = L & ~same
        res = bool(L.any()); self.cache[key] = res; return res
    def simple(self, ci): return not self.crosses(ci, ci)
    def parity(self, ci, x, y, depth):
        sel = self.lens <= comp_len(self.comps[ci]) + depth + 2
        A, B = self.lifts(ci)
        return int(linked_vec(x, y, A[sel], B[sel]).sum()) % 2

def run_config(C, rng, eps_list=(1e-3, 1e-5), verbose=False):
    comps = build_components(C)
    topo = Topology(C, comps, 0)
    good = [ci for ci in range(len(comps)) if topo.simple(ci)]
    n = len(C.marks); N = C.N
    # sample points: one per (arc,piece segment); and base point
    samples = []
    for pid, segs in enumerate(C.pieces):
        for (k, f0, f1) in segs:
            if f1 - f0 > 1e-9:
                samples.append((arc_point(C.arcs, k, (f0 + f1)/2), pid))
    s0, p0 = samples[0]
    S1 = C.mets[0]
    gen_imgs = [Surface.act(S1.gens[k], s0) for k in (0, 2)]
    def region_bits(m):
        return [(C.labels[pid] < N and (m >> C.labels[pid]) & 1) for pid in range(len(C.pieces))]
    def mark_endpoints(m):
        bits = region_bits(m)
        return [mi for mi in range(n) if bits[C.adj_pieces[mi][0]] != bits[C.adj_pieces[mi][1]]]
    arcs_at = {mi: [ci for ci in good if comps[ci][0] == 'arc' and comps[ci][1] != comps[ci][3] and mi in (comps[ci][1], comps[ci][3])] for mi in range(n)}
    closed = [ci for ci in good if comps[ci][0] == 'closed']
    K = len(C.mets)
    w = C.w
    def cost(ci, eps, mode='sum'):
        ls = [C.length(comps[ci], j, eps) for j in range(K)]
        return ls
    results = {}
    for eps in eps_list:
        Lc = {ci: np.array(cost(ci, eps)) for ci in good}
        F = {}; Fmax = {}; arg = {}
        for m in range(1, 2**N):
            E = mark_endpoints(m)
            bits = region_bits(m)
            best = [np.inf, None]; bestmax = [np.inf]
            # enumerate matchings with compatible arcs
            def rec(unmatched, chosen, acc):
                if not unmatched:
                    # closed subsets (size <= 3), check homology
                    cl_ok = [ci for ci in closed if all(not topo.crosses(ci, a) for a in chosen)]
                    for r in range(0, 4):
                        for sub in itertools.combinations(cl_ok, r):
                            if any(topo.crosses(x, y) for x, y in itertools.combinations(sub, 2)): continue
                            cut = list(chosen) + list(sub)
                            # homology parity
                            ok = True
                            for (s, pid) in samples:
                                par = sum(topo.parity(ci, s0, s, 0) for ci in cut) % 2
                                if par != (int(bits[pid]) ^ int(bits[p0])): ok = False; break
                            if ok:
                                for y in gen_imgs:
                                    if sum(topo.parity(ci, s0, y, 1) for ci in cut) % 2: ok = False; break
                            if not ok: continue
                            tot = sum(Lc[ci] for ci in cut) if cut else np.zeros(K)
                            v = float(np.dot(w, tot))
                            if v < best[0] - 1e-12: best[0], best[1] = v, tuple(cut)
                            bestmax[0] = min(bestmax[0], float(np.max(tot)) if cut else 0.0)
                    return
                p = unmatched[0]
                for ci in arcs_at[p]:
                    q = comps[ci][3] if comps[ci][1] == p else comps[ci][1]
                    if q not in unmatched[1:]: continue
                    if any(topo.crosses(ci, a) for a in chosen): continue
                    rec([x for x in unmatched[1:] if x != q], chosen + [ci], acc)
            rec(E, [], 0.0)
            F[m] = best[0]; Fmax[m] = bestmax[0]; arg[m] = best[1]
        # t_p
        t = np.zeros(n)
        for mi in range(n):
            for j in range(K):
                Sj = C.mets[j]; p = C.pos[j][mi]
                imgs = []
                for wd, M in zip(C.W_CHK, C.mats[j]):
                    if len(wd) > 5: continue
                    for k in range(n):
                        z = Surface.act(M, C.pos[j][k])
                        if abs(z - p) > 1e-12: imgs.append(-1.0/(z - p))
                d = max(imgs) - min(imgs)
                t[mi] += w[j]*np.log(1.0/(eps*C.hscale[j][mi]*d))
        G = {m: F[m] - sum(t[mi] for mi in mark_endpoints(m)) for m in F}
        results[eps] = dict(F=F, G=G, Fmax=Fmax, t=t, arg=arg)
    return results, comps, good

def main(topo, N, trials, seed, kmet=2):
    rng = np.random.default_rng(seed)
    done = 0; stats = dict(dG=0.0, fF=np.inf, fG=np.inf, fGctl=np.inf, fmax=np.inf, nviol_ctl=0, nviol_max=0,
                           inf=0, tmin=np.inf, stab=0.0)
    names = {}
    while done < trials:
        nmark = int(rng.integers(2, 6)) if topo == 'pants' else int(rng.integers(2, 6))
        nmark += nmark % 2 if topo == 'torus' else 0
        C = Config(rng, topo, kmet, N, nmark)
        C.w = rng.uniform(0.2, 1.0, size=kmet)
        if not C.assign(rng): continue
        res, comps, good = run_config(C, rng)
        r1, r2 = res[1e-3], res[1e-5]
        if any(not np.isfinite(v) for v in r1['F'].values()):
            stats['inf'] += 1; continue
        done += 1
        stats['dG'] = max(stats['dG'], max(abs(r1['G'][m] - r2['G'][m]) for m in r1['G']))
        vF, nF = fc.facet_check(r2['F'], N); vG, nG = fc.facet_check(r1['G'], N)
        stats['fF'] = min(stats['fF'], vF); stats['fG'] = min(stats['fG'], vG)
        if vG < -1e-9: names[nG] = names.get(nG, 0) + 1
        stats['tmin'] = min(stats['tmin'], float(r2['t'].min()))
        # controls
        W = C.w.sum()
        bp = lambda m: 0
        Gc = {m: r1['G'][m] - 0.5*W*(r1['F'][m] - r1['G'][m])/max(1e-300, 1) * 0 for m in r1['G']}
        # shifted control: G - 0.5*W*sum_{p in dA} 1
        n = len(C.marks)
        def nends(m):
            bits = [(C.labels[pid] < N and (m >> C.labels[pid]) & 1) for pid in range(len(C.pieces))]
            return sum(bits[C.adj_pieces[mi][0]] != bits[C.adj_pieces[mi][1]] for mi in range(n))
        Gc = {m: r1['G'][m] - 0.5*W*nends(m) for m in r1['G']}
        vc, _ = fc.facet_check(Gc, N); stats['fGctl'] = min(stats['fGctl'], vc); stats['nviol_ctl'] += vc < -1e-9
        Fm = {m: r2['Fmax'][m] - W*0 for m in r2['Fmax']}
        # max functional with endpoint divergences removed the same way (use t of a single metric scale)
        vm, _ = fc.facet_check({m: r2['Fmax'][m] - (r2['F'][m] - r2['G'][m])/W for m in r2['Fmax']}, N)
        stats['fmax'] = min(stats['fmax'], vm); stats['nviol_max'] += vm < -1e-9
        print(f"[{topo} N={N}] {done}/{trials} marks={n} pieces={len(C.pieces)} comps={len(good)} "
              f"dG={stats['dG']:.1e} facet(F)={vF:+.1e} facet(G)={vG:+.1e} ctl={vc:+.1e} max={vm:+.1e}", flush=True)
    print("SUMMARY", topo, N, stats, names, flush=True)

if __name__ == "__main__":
    topo = sys.argv[1]; N = int(sys.argv[2]); trials = int(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    main(topo, N, trials, seed)
