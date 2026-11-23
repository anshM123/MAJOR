r"""
rh_lemma_count (R4ARC_HEFT.md sec. 2.3, Lemma Y): exhaustive check of the counting lemma for spinning contact diagrams of four spin-2
operators (Delta = 3, d = 3) written as  (monomial in embedding invariants [x one epsilon]) x D_{Delta'}(P),  Delta'_i = 3 + deg_{P_i}.
In the og light-ray configuration (x1 = -t1 z1, x3 = t3 z3, Z1 = z1, Z3 = z3) the invariants Z1.Z3, Z1.P3, Z3.P1, P1.P3 are O(y), all
others O(1); an epsilon with >= 3 of {P1, Z1, P3, Z3} is O(y^{1/2}) (exactly 0 with all four).  The D-function's leftmost
(13)-channel Mellin pole sits at d13 = max(0, sigma), sigma = (D'1 + D'3 - D'2 - D'4)/2, i.e. y^{-max(0, sigma)}.
Claim: m := (y-order of the prefactor) >= max(0, sigma) for every such term.  Enumerated: every assignment of the 8 Z-legs (two per
Z_i) to targets, P.P factors up to total multiplicity 2, and parity-even (no epsilon) or one epsilon on any 5 of the 8 vectors.
Also verifies the y-orders numerically (finite differences in log y) for a random sample of the parity-odd epsilons.
"""
import itertools
import numpy as np

V = ['P1', 'P2', 'P3', 'P4', 'Z1', 'Z2', 'Z3', 'Z4']
FAM13 = {'P1', 'Z1', 'P3', 'Z3'}


def is13(u, v):
    return {u, v} in ({'Z1', 'Z3'}, {'Z1', 'P3'}, {'Z3', 'P1'}, {'P1', 'P3'})


def check():
    worst = -99; nterms = 0
    # targets for a Z_i leg: other Z_j or P_j (j != i)
    Zs = ['Z1', 'Z2', 'Z3', 'Z4']
    for eps in [None] + list(itertools.combinations(V, 5)):
        legs_used = {z: (1 if eps and z in eps else 0) for z in Zs}
        if any(v > 2 for v in legs_used.values()):
            continue
        free = []
        for z in Zs:
            free += [z]*(2 - legs_used[z])
        # pair Z-legs: each free leg goes to a P_j (j != i) or is paired with another free Z-leg of a different Z
        def assign(legs):
            if not legs:
                yield []
                return
            a = legs[0]; rest = legs[1:]
            i = a[1]
            for Pj in ['P1', 'P2', 'P3', 'P4']:
                if Pj[1] != i:
                    for r in assign(rest):
                        yield [(a, Pj)] + r
            for k, b in enumerate(rest):
                if b != a:
                    for r in assign(rest[:k] + rest[k + 1:]):
                        yield [(a, b)] + r
        seen = set()
        for pairs in assign(free):
            key = tuple(sorted(tuple(sorted(p)) for p in pairs))
            if key in seen:
                continue
            seen.add(key)
            for PP in [()] + [(p,) for p in itertools.combinations(['P1', 'P2', 'P3', 'P4'], 2)] + \
                      [pp for pp in itertools.combinations_with_replacement(list(itertools.combinations(['P1', 'P2', 'P3', 'P4'], 2)), 2)]:
                PP = [p for p in PP] if PP and isinstance(PP[0], tuple) else []
                deg = {p: 0 for p in ['P1', 'P2', 'P3', 'P4']}
                m = 0.0
                for (u, v) in pairs:
                    for w in (u, v):
                        if w in deg:
                            deg[w] += 1
                    m += is13(u, v)
                for (u, v) in PP:
                    deg[u] += 1; deg[v] += 1
                    m += is13(u, v)
                if eps:
                    for w in eps:
                        if w in deg:
                            deg[w] += 1
                    k13 = len(FAM13 & set(eps))
                    if k13 == 4:
                        continue                      # identically zero
                    if k13 == 3:
                        m += 0.5
                sigma = (deg['P1'] + deg['P3'] - deg['P2'] - deg['P4'])/2
                worst = max(worst, max(0, sigma) - m)
                nterms += 1
    return nterms, worst


def numeric_eps(ntry=200, seed=3):
    """y-order of epsilon(5 of 8) in the light-ray configuration, measured numerically."""
    rng = np.random.default_rng(seed)
    def emb(x, z=None):
        x = np.asarray(x, float)
        if z is None:
            return np.array([1.0, x @ np.diag([-1, 1, 1]) @ x, *x])
        z = np.asarray(z, float)
        return np.array([0.0, 2*(x @ np.diag([-1, 1, 1]) @ z), *z])
    bad = 0; out = []
    for _ in range(ntry):
        t1, t3 = rng.uniform(0.2, 2, 2); phib = rng.uniform(0, 6.28)
        x2 = rng.normal(size=3); x4 = rng.normal(size=3)
        def nullv():
            v = rng.normal(size=2); v /= np.linalg.norm(v); return np.array([1.0, *v])
        z2, z4 = nullv(), nullv()
        sub = list(itertools.combinations(range(8), 5))[rng.integers(0, 56)]
        vals = []
        for dphi in (1e-2, 1e-3):
            z1 = np.array([1, np.cos(phib + dphi/2), np.sin(phib + dphi/2)]); z3 = np.array([1, np.cos(phib - dphi/2), np.sin(phib - dphi/2)])
            vecs = [emb(-t1*z1), emb(x2), emb(t3*z3), emb(x4), emb(-t1*z1, z1), emb(x2, z2), emb(t3*z3, z3), emb(x4, z4)]
            M = np.array([vecs[i] for i in sub])
            vals.append(abs(np.linalg.det(M))/np.prod(np.linalg.norm(M, axis=1)))
        y1, y2 = 2*np.sin(0.5e-2)**2, 2*np.sin(0.5e-3)**2
        k13 = len(FAM13 & {V[i] for i in sub})
        if max(vals) < 1e-12:
            order = np.inf                                        # identically zero (roundoff)
        else:
            order = np.log(vals[0]/vals[1])/np.log(y1/y2)
        # claim: order >= 0 for k13 <= 2 (bounded), = 1/2 for k13 = 3, identically 0 for k13 = 4
        if (k13 <= 2 and order < -0.2) or (k13 == 3 and abs(order - 0.5) > 0.05) or (k13 == 4 and np.isfinite(order)):
            bad += 1
        out.append((k13, order))
    return bad, out


if __name__ == '__main__':
    n, w = check()
    lines = [f'Lemma Y exhaustive count: {n} (monomial x [epsilon]) structures; max over all of [max(0, sigma) - m] = {w:+.2f} (claim: <= 0)']
    bad, out = numeric_eps()
    by = {}
    for k13, o in out:
        by.setdefault(k13, []).append(o)
    lines.append('numerical y-order of epsilon(5 vectors) in the light-ray configuration, by number k13 of {P1,Z1,P3,Z3} among them:')
    for k13 in sorted(by):
        arr = np.array(by[k13]); fin = arr[np.isfinite(arr)]
        lines.append(f'   k13 = {k13}: n = {len(arr)}, measured order min/max = {fin.min() if fin.size else np.nan:.3f} / {fin.max() if fin.size else np.nan:.3f}'
                     f' (claim: {">= 0 up to O(dphi) drift (bounded)" if k13 < 3 else ("1/2" if k13 == 3 else "identically 0")}); non-finite (exact zero): {np.sum(~np.isfinite(arr))}')
    lines.append(f'   mismatches: {bad}')
    print('\n'.join(lines))
    open('cache/lemmaY.txt', 'w').write('\n'.join(lines) + '\n')
