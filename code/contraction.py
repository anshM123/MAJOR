"""Holographic entropy inequalities and contraction maps.

An inequality is (lhs, rhs): lists of (coefficient, frozenset_of_parties), meaning
    sum_l alpha_l S(A_l) >= sum_r beta_r S(B_r).
find_contraction_map searches for f:{0,1}^L -> {0,1}^R with
    ||f(x)-f(x')||_beta <= ||x-x'||_alpha   and   f(x^(i)) = y^(i) for all parties i and the purifier.
A returned map is an exact certificate (integer arithmetic) and is re-verified by verify()."""
import itertools

def occurrence(terms, i): return tuple(1 if i in A else 0 for _, A in terms)

def find_contraction_map(lhs, rhs, parties):
    L, R = len(lhs), len(rhs); al = [c for c, _ in lhs]; be = [c for c, _ in rhs]
    X = list(itertools.product((0, 1), repeat=L)); Y = list(itertools.product((0, 1), repeat=R))
    dX = lambda x, y: sum(a for a, s, t in zip(al, x, y) if s != t)
    dY = lambda x, y: sum(b for b, s, t in zip(be, x, y) if s != t)
    fixed = {occurrence(lhs, i): occurrence(rhs, i) for i in list(parties) + [None]}
    for a, b in itertools.combinations(fixed, 2):
        if dY(fixed[a], fixed[b]) > dX(a, b): return None
    order = sorted(X, key=lambda x: (x not in fixed, min(dX(x, a) for a in fixed)))
    f = {}
    def cands(x):
        if x in fixed: return [fixed[x]]
        return [y for y in Y if all(dY(y, f[z]) <= dX(x, z) for z in f)]
    def rec(k):
        if k == len(order): return True
        x = order[k]
        for y in cands(x):
            if all(dY(y, f[z]) <= dX(x, z) for z in f):
                f[x] = y
                if rec(k+1): return True
                del f[x]
        return False
    return dict(f) if rec(0) else None

def verify(lhs, rhs, parties, f):
    al = [c for c, _ in lhs]; be = [c for c, _ in rhs]
    for x in f:
        for z in f:
            if sum(b for b, s, t in zip(be, f[x], f[z]) if s != t) > sum(a for a, s, t in zip(al, x, z) if s != t): return False
    return all(f[occurrence(lhs, i)] == occurrence(rhs, i) for i in list(parties) + [None]) and len(f) == 2**len(lhs)

def S(*ps): return frozenset(ps)
def cyclic(k):
    """sum_i S(A_i..A_{i+k}) >= sum_i S(A_{i+1}..A_{i+k}) + S(A_1..A_{2k+1})   (k=1: MMI)"""
    n = 2*k+1
    lhs = [(1, frozenset((i+j) % n for j in range(k+1))) for i in range(n)]
    rhs = [(1, frozenset((i+1+j) % n for j in range(k))) for i in range(n)] + [(1, frozenset(range(n)))]
    return lhs, rhs, list(range(n))
SSA = ([(1, S(0, 1)), (1, S(1, 2))], [(1, S(1)), (1, S(0, 1, 2))], [0, 1, 2])
# a deliberately FALSE "inequality" used as a control: reversed MMI
REV_MMI = ([(1, S(0)), (1, S(1)), (1, S(2)), (1, S(0, 1, 2))], [(1, S(0, 1)), (1, S(1, 2)), (1, S(0, 2))], [0, 1, 2])

def catalogue():
    out = {"SSA": SSA, "MMI": cyclic(1), "cyclic5": cyclic(2), "cyclic7": cyclic(3)}
    return out

if __name__ == "__main__":
    for name, (lhs, rhs, parties) in list(catalogue().items()) + [("reversed-MMI (control)", REV_MMI)]:
        f = find_contraction_map(lhs, rhs, parties)
        print(f"{name:24s} L={len(lhs)} R={len(rhs)}: contraction map", "FOUND and verified" if f and verify(lhs, rhs, parties, f) else "does not exist / not found")
