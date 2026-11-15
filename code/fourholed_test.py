"""Falsification test of Theorem 1 in a multiply connected topology: the four-holed sphere (pi_1 = F_3),
regions = unions of whole boundary components (so the inequality is exact, no cut-off).

Hyperbolic structures: classical Schottky groups <g1,g2,g3>, g_i pairing adjacent disjoint intervals I_i, I_i' of the
real line (order I_1,I_1',I_2,I_2',I_3,I_3').  Boundary classes: a=g3, b=g2, c=g1, d=(g3 g2 g1)^{-1}.
Essential simple closed curves: x1x2 and x2x3 for all images (x1,x2,x3) of (a,b,c) under the braid group B_3
(sigma_1: (x1,x2,x3)->(x1 x2 x1^-1, x1, x3), sigma_2 similarly), enumerated to a given depth.
Cut classes for a set T of boundaries: a set B of boundary geodesics, optionally together with ONE essential curve e
(two distinct essential curves on the four-holed sphere always intersect), such that B (symmetric difference) side(e)
equals T or its complement.   F(T) = min over classes of sum_j w_j * length_j(class)   -- same class in all metrics.
k=2 with w=(1/2,1/2) is the HRT entropy of the rotating four-boundary vacuum wormhole."""
import numpy as np, sys, itertools

def mob(l, r):  # affine map taking [l,r] to [-1,1], as a matrix
    return np.array([[2.0, -(l+r)], [0.0, r-l]])
def schottky(rng):
    # 6 disjoint intervals; wide intervals with narrow gaps give short boundary/neck geodesics
    widths = rng.uniform(0.3, 3.0, size=6); gaps = np.exp(rng.uniform(np.log(0.02), np.log(2.0), size=6))
    pos = 0.0; I = []
    for wd, gp in zip(widths, gaps):
        I.append((pos, pos+wd)); pos += wd+gp
    gens = []
    for i in range(3):
        (l, r), (lp, rp) = I[2*i], I[2*i+1]; t = rng.normal()*0.7
        M = np.array([[np.cosh(t), np.sinh(t)], [np.sinh(t), np.cosh(t)]]); S = np.array([[0.0, -1.0], [1.0, 0.0]])
        g = np.linalg.inv(mob(lp, rp)) @ M @ S @ mob(l, r)
        g = g/np.sqrt(np.linalg.det(g)); gens.append(g)
    return gens   # g1,g2,g3

def pants(l1, l2, l3):
    """A, B in SL(2,R) with tr A = -2cosh(l1/2), tr B = -2cosh(l2/2), tr AB = -2cosh(l3/2): a pair-of-pants group"""
    lam = np.exp(l1/2); c2 = np.cosh(l2/2); c3 = np.cosh(l3/2)
    A = -np.diag([lam, 1/lam]); p = (2*c3 + 2*c2/lam)/(lam - 1/lam); sB = -2*c2 - p
    B = np.array([[p, 1.0], [p*sB - 1.0, sB]]); return A, B
def fenchel_nielsen(la, lb, lc, ld, le, tw):
    """four-holed sphere: pants (a,b,e) and (c,d,e^-1) glued along the neck e of length le with twist tw; abcd=1"""
    A, B = pants(la, lb, le); E = A @ B; Cp, Dp = pants(lc, ld, le)
    Rm = np.diag([1.0, -1.0]); Cp, Dp = Rm @ Cp @ Rm, Rm @ Dp @ Rm   # mirror image: puts the second pair of pants on the
    # other side of the neck axis (validated numerically: with this choice every enumerated simple curve is hyperbolic,
    # without it elliptic elements appear)
    def eig(M):
        w, V = np.linalg.eig(M); o = np.argsort(-np.abs(w.real)); return V[:, o].real
    U = eig(inv(E)); V = eig(Cp @ Dp)
    g = U @ np.diag([np.exp(tw/2), np.exp(-tw/2)]) @ np.linalg.inv(V); g = g/np.sqrt(abs(np.linalg.det(g)))
    gi = inv(g)*np.sign(np.linalg.det(g))
    C = g @ Cp @ gi; D = g @ Dp @ gi
    res = np.abs(A @ B @ C @ D - np.eye(2)).max()
    if res > 1e-9: return None          # badly conditioned sample: reject (caller resamples)
    return A, B, C
def sample_fn(rng):
    la, lb, lc, ld = rng.uniform(1.0, 7.0, size=4); le = np.exp(rng.uniform(np.log(0.3), np.log(5.0)))
    while True:
        la, lb, lc, ld = rng.uniform(1.0, 6.0, size=4); le = np.exp(rng.uniform(np.log(0.3), np.log(4.0)))
        tw = rng.uniform(-2.0, 2.0)*le
        out = fenchel_nielsen(la, lb, lc, ld, le, tw)
        if out is not None: return out

def length(M):
    with np.errstate(all='ignore'):
        t = np.abs(np.trace(M, axis1=-2, axis2=-1)); L = 2*np.arccosh(np.maximum(t/2, 1.0))
    return np.where(np.isfinite(L), L, np.inf)   # overflowed (astronomically long) words can never be minimisers
def inv(M):
    """exact inverse for determinant-one 2x2 matrices (stable for long words)"""
    R = np.empty_like(M); R[..., 0, 0] = M[..., 1, 1]; R[..., 1, 1] = M[..., 0, 0]; R[..., 0, 1] = -M[..., 0, 1]; R[..., 1, 0] = -M[..., 1, 0]
    return R

EPS = 2.3e-16
def _norm(M): return np.sqrt((M*M).sum(axis=(-1, -2)))
def _mul(A, EA, B, EB):
    """product with a running bound on the absolute rounding error of the entries"""
    C = A @ B; nA = _norm(A); nB = _norm(B)
    return C, EA*nB + nA*EB + 4*EPS*nA*nB

def enumerate_curves(triples, depth, tol=1e-7):
    """triples: array (k,3,2,2) for k metrics.  Returns pair types (n,2) and lengths (n,k).
    Every matrix carries a rigorous-style bound on its accumulated rounding error; a curve whose trace is not
    known to absolute accuracy `tol` in some metric is discarded (length = +inf).  Such curves are astronomically
    long words; the driver checks that no minimiser sits at the boundary of the reliable set."""
    np.seterr(all='ignore')
    X = triples[None].copy(); E = np.zeros(X.shape[:3])
    types = np.array([[0, 1, 2]]); last = np.array([-1]); out_t = []; out_l = []; out_d = []
    def record(X, E, types, d):
        for (i, j) in ((0, 1), (1, 2)):
            P, EP = _mul(X[:, :, i], E[:, :, i], X[:, :, j], E[:, :, j])
            L = length(P); L = np.where(2*EP < tol, L, np.inf)
            out_l.append(L); out_t.append(np.sort(types[:, [i, j]], axis=1)); out_d.append(np.full(len(X), d))
    record(X, E, types, 0)
    for d in range(1, depth+1):
        newX = []; newE = []; newT = []; newL = []
        for mv in range(4):
            keep = last != (mv ^ 1)
            if not keep.any(): continue
            Y = X[keep].copy(); EY = E[keep].copy(); T = types[keep].copy(); i = 0 if mv < 2 else 1
            xi = Y[:, :, i].copy(); xj = Y[:, :, i+1].copy(); ei = EY[:, :, i].copy(); ej = EY[:, :, i+1].copy()
            if mv % 2 == 0:
                P, EP = _mul(xi, ei, xj, ej); P, EP = _mul(P, EP, inv(xi), ei)
                Y[:, :, i] = P; EY[:, :, i] = EP; Y[:, :, i+1] = xi; EY[:, :, i+1] = ei
            else:
                P, EP = _mul(inv(xj), ej, xi, ei); P, EP = _mul(P, EP, xj, ej)
                Y[:, :, i] = xj; EY[:, :, i] = ej; Y[:, :, i+1] = P; EY[:, :, i+1] = EP
            T[:, [i, i+1]] = T[:, [i+1, i]]
            newX.append(Y); newE.append(EY); newT.append(T); newL.append(np.full(len(Y), mv))
        X = np.concatenate(newX); E = np.concatenate(newE); types = np.concatenate(newT); last = np.concatenate(newL)
        record(X, E, types, d)
    return np.concatenate(out_t), np.concatenate(out_l), np.concatenate(out_d)

def entropies(metrics, w, depth, mode="sum"):
    """metrics: list of (a,b,c) matrix triples.  returns dict T -> F(T) for T subset of {0,1,2} (3 = purifier d)"""
    k = len(metrics); tr = np.array([[m[0], m[1], m[2]] for m in metrics])          # (k,3,2,2)
    bl = np.array([[length(m[0]), length(m[1]), length(m[2]), length(inv(m[0] @ m[1] @ m[2]))] for m in metrics])  # (k,4)
    types, L, dep = enumerate_curves(tr, depth)
    comb = (L*w).sum(1) if mode == "sum" else L.max(1)
    comb = np.where(np.isfinite(comb), comb, np.inf)
    best = {}; arg = {}
    for pair in ((0, 1), (0, 2), (1, 2)):
        sel = (types[:, 0] == pair[0]) & (types[:, 1] == pair[1])
        idx = np.where(sel)[0]; j = idx[np.argmin(comb[sel])]; best[pair] = L[j]; arg[pair] = j
    def cost(B, e=None):
        v = bl[:, list(B)].sum(1) if B else np.zeros(k)
        if e is not None: v = v + best[e]
        return (v*w).sum() if mode == "sum" else v.max()
    full = frozenset(range(4)); F = {}; used_essential = False
    for r in range(0, 4):
        for T in itertools.combinations(range(3), r):
            T = frozenset(T); c0 = min(cost(T), cost(full - T)); c = c0
            for pair in best:
                side = frozenset(pair)
                c = min(c, cost(T ^ side, pair), cost((full - T) ^ side, pair))
            F[T] = c; used_essential = used_essential or (c < c0 - 1e-12)
    return F, {p: int(dep[j]) for p, j in arg.items() if np.isfinite(comb[j])}, used_essential

SAMPLER = "fn"
if __name__ == "__main__":
    if len(sys.argv) > 4: SAMPLER = sys.argv[4]
    seed = int(sys.argv[1]); trials = int(sys.argv[2]); depth = int(sys.argv[3]); rng = np.random.default_rng(seed)
    fs = lambda *x: frozenset(x)
    for k, mode in ((1, "sum"), (2, "sum"), (3, "sum"), (2, "max")):
        worst = np.inf; viol = 0; deep = 0; fracs = []
        for t in range(trials):
            metrics = []
            for j in range(k):
                metrics.append(sample_fn(rng) if SAMPLER == 'fn' else (lambda g: (g[2], g[1], g[0]))(schottky(rng)))
            w = np.full(k, 1.0/k) if k == 2 else rng.uniform(0.2, 1, size=k)
            F, argdepth, frac = entropies(metrics, w, depth, mode)
            mmi = F[fs(0, 1)] + F[fs(1, 2)] + F[fs(0, 2)] - F[fs(0)] - F[fs(1)] - F[fs(2)] - F[fs(0, 1, 2)]
            ssa = min(F[fs(0, 1)] + F[fs(1, 2)] - F[fs(1)] - F[fs(0, 1, 2)], F[fs(0, 1)] + F[fs(0, 2)] - F[fs(0)] - F[fs(0, 1, 2)])
            sc = sum(abs(v) for v in F.values()) + 1
            if min(mmi, ssa) < -1e-9*sc: viol += 1
            worst = min(worst, mmi/sc, ssa/sc)
            # was any minimising curve found only at the last BFS level? (then depth may be insufficient)
            deep = max(deep, max(list(argdepth.values()) + [0])); fracs.append(frac)
        print(f"k={k} mode={mode}: trials={trials} depth={depth}: violations of MMI/SSA = {viol}; min normalised (LHS-RHS) = {worst:+.2e}; deepest braid level of any minimiser = {deep}; fraction of trials where an essential curve is part of some minimal cut = {np.mean(fracs):.2f}", flush=True)
