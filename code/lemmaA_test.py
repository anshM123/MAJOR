"""Direct falsification test of Lemma A (the topological contraction lemma) on the cylinder with marked points.

Everything is exact integer combinatorics.  Marked points 0,1,...,n-1 on one ideal circle (period n); elementary arc k
is (k,k+1); the other end of the cylinder is a purifier piece.  Universal cover of the marked circle = the real line.
 cut components:  chord i<j, hole outside  -> lift (i, j);   hole inside -> lift (j, i+n);   core curve.
 test classes:    T1(u1,u1+m): arc from piece u1 to the piece at u1+m (m>=1, may wind several times);  T2(u): arc
                  from piece u to the other end;  intersection numbers by linking of lifts:
   i(T1(u1,u2), chord(t1,t2)) = #{k : exactly one of u1+kn, u2+kn lies in (t1,t2)},  i(T1, core) = 0,
   i(T2(u),   chord(t1,t2)) = #{k : u+kn in (t1,t2)},                                 i(T2, core) = 1.
Test (L1): for every single test class tau, F_tau(A) = min_c i(tau,c) satisfies each inequality  (consequence of A1).
Test (L2): for ARBITRARY classes gamma_l in C(A_l) -- not minimisers of anything -- there exist delta_r in C(B_r) with
           sum_r beta_r i(tau,delta_r) <= sum_l alpha_l i(tau,gamma_l)  for ALL test classes tau simultaneously
           (exhaustive search over tuples of classes).  This is statement (A1) itself."""
import numpy as np, sys, itertools
from annulus_test import region_classes
from contraction import catalogue, REV_MMI

def tests(n, maxwind):
    T = [("T1", u+0.5, u+0.5+m) for u in range(n) for m in range(1, maxwind*n+1)] + [("T2", u+0.5, None) for u in range(n)]
    return T
def inter_chord(test, t1, t2, n):
    kind, u1, u2 = test; cnt = 0
    for k in range(-8, 9):
        a = t1 < u1 + k*n < t2
        if kind == "T2": cnt += a
        else: cnt += (a != (t1 < u2 + k*n < t2))
    return cnt
def class_vectors(Aarcs, n, T):
    out = []
    for M, sig, core in region_classes(frozenset(Aarcs), n):
        v = np.zeros(len(T), int)
        for (i, j), inside in zip(M, sig):
            t1, t2 = (j, i+n) if inside else (i, j)
            v += np.array([inter_chord(t, t1, t2, n) for t in T])
        if core: v += np.array([1 if t[0] == "T2" else 0 for t in T])
        out.append(v)
    return np.array(out)                      # (#classes, #tests)

def exists_dominating(rhs_lists, betas, budget):
    """DFS: choose one row from each list so that sum_r beta_r * row_r <= budget componentwise"""
    order = sorted(range(len(rhs_lists)), key=lambda r: len(rhs_lists[r]))
    def rec(idx, rem):
        if idx == len(order): return True
        r = order[idx]
        for row in rhs_lists[r]:
            nr = rem - betas[r]*row
            if (nr >= 0).all() and rec(idx+1, nr): return True
        return False
    return rec(0, budget.copy())

if __name__ == "__main__":
    seed = int(sys.argv[1]); trials = int(sys.argv[2]); rng = np.random.default_rng(seed); cat = catalogue()
    for name in ["SSA", "MMI", "cyclic5", "cyclic7", "revMMI(control)"]:
        lhs, rhs, parties = REV_MMI if name.startswith("rev") else cat[name]; N = len(parties)
        L1_viol = 0; L1_checked = 0; L2_fail = 0; L2_checked = 0
        for t in range(trials if N <= 5 else trials//4):
            n = int(rng.integers(N+1, N+4))
            while True:
                lab = rng.integers(-1, N, size=n)
                if set(range(N)) <= set(lab.tolist()) and -1 in lab: break
            T = tests(n, 2)
            vecs = {}
            def cv(A):
                A = frozenset(A)
                if A not in vecs: vecs[A] = class_vectors([a for a in range(n) if lab[a] in A], n, T)
                return vecs[A]
            Lv = [cv(A) for _, A in lhs]; Rv = [cv(A) for _, A in rhs]
            al = np.array([c for c, _ in lhs]); be = np.array([c for c, _ in rhs])
            # (L1) single test classes
            FL = np.array([v.min(axis=0) for v in Lv]); FR = np.array([v.min(axis=0) for v in Rv])
            q = al @ FL - be @ FR; L1_viol += int((q < 0).sum()); L1_checked += q.size
            # (L2) arbitrary gamma tuples
            for rep in range(6):
                budget = sum(a*v[rng.integers(len(v))] for a, v in zip(al, Lv))
                ok = exists_dominating(Rv, be, budget); L2_fail += (not ok); L2_checked += 1
        print(f"{name:16s}: (L1) single-test-class violations {L1_viol}/{L1_checked};   (L2) gamma-tuples with NO dominating delta-tuple {L2_fail}/{L2_checked}", flush=True)
