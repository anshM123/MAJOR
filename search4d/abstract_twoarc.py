"""Abstract tests with the copies engine.

 T0  copies engine == one-arc DP engine (abstract_axis.entropy_vector).
 T1  two free arcs (two-sided SO(2)-symmetric geometry: bdry1, S-axis, bdry2,
     N-axis): circular-decomposable augmented metric -> facets?
 T2  two free arcs, generic (non-Kalmanson) perturbations of T1 data, kept only
     if ALL physical SSA and weak-monotonicity instances hold (all unions of
     asymptotic arcs, pure total state): any facet violation?
Run: python3 -I abstract_twoarc.py [nsamp]
"""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from abstract_axis import F, ORB, SUBSETS, random_labels, entropy_vector, circ_decomp
from circle_engine import Circle, INF


def onearc_seq(labels, m):
    seq = []
    for j in range(m):
        seq += [('arc', labels[j]), ('edge', j)]
    seq += [('arc', labels[m]), ('free', 0)]
    return seq


def twoarc_seq(lab1, lab2):
    """lab1: labels of bdry1 arcs J_0..J_p (N->S), lab2: K_0..K_q (S->N)."""
    p = len(lab1) - 1; q = len(lab2) - 1
    seq = []
    for j in range(p):
        seq += [('arc', lab1[j]), ('edge', j)]
    seq += [('arc', lab1[p]), ('free', 1)]          # S axis = free arc 1
    for j in range(q):
        seq += [('arc', lab2[j]), ('edge', p + j)]
    seq += [('arc', lab2[q]), ('free', 0)]          # N axis = free arc 0
    return seq


def twoarc_circ(p, q, rng, sparse=0.3):
    # augmented circle order: N, a_1..a_p, S, b_1..b_q
    n = p + q + 2
    order = ['N'] + [('e', j) for j in range(p)] + ['S'] + [('e', p + j) for j in range(q)]
    D = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            if rng.random() < sparse:
                ins = np.zeros(n, bool); ins[a:b] = True
                D += rng.exponential() * (ins[:, None] != ins[None, :])
    iN, iS = 0, p + 1
    eidx = [k for k in range(n) if k not in (iN, iS)]
    L = D[np.ix_(eidx, eidx)]
    C = [D[iN, eidx], D[iS, eidx]]
    H = np.array([[0, D[iN, iS]], [D[iN, iS], 0]])
    return L, C, H


def vec_from(circ):
    return np.array([circ.S(set(S)) for S in SUBSETS])


def physical_slack(seq_builder, nar, L, C, H, k):
    """min SSA / WM slack over all pairs of unions of asymptotic arcs (each arc its own label)."""
    circ = Circle(seq_builder(list(range(nar))), L, C, H, k=k)
    allarcs = frozenset(range(nar))
    S = {}
    for r in range(0, nar + 1):
        for U in itertools.combinations(range(nar), r):
            U = frozenset(U)
            S[U] = circ.S(U) if U else 0.0
    keys = list(S)
    worst = np.inf
    for U in keys:
        for V in keys:
            ssa = S[U] + S[V] - S[U | V] - S[U & V]
            wm = S[U] + S[V] - S[U - V] - S[V - U]
            worst = min(worst, ssa, wm)
    # purity check (should be 0)
    pur = max(abs(S[U] - S[allarcs - U]) for U in keys)
    return worst, pur


if __name__ == '__main__':
    nsamp = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    rng = np.random.default_rng(7)
    # T0
    dev = 0
    for s in range(20):
        m = int(rng.integers(5, 8)); L, c = circ_decomp(m, rng)
        lab = random_labels(m, rng)
        v1 = entropy_vector(lab, L, c)
        circ = Circle(onearc_seq(lab, m), L, [c], np.zeros((1, 1)), k=4)
        dev = max(dev, np.abs(vec_from(circ) - v1).max())
    print('T0 max |copies - DP| =', dev, flush=True)
    # T1
    worst = np.full(len(ORB), np.inf)
    for s in range(nsamp):
        p = int(rng.integers(2, 5)); q = int(rng.integers(2, 5))
        L, C, H = twoarc_circ(p, q, rng)
        for _ in range(5):
            lab = [int(x) for x in rng.integers(0, 6, size=p + q + 2)]
            if any(lab[j] == lab[j + 1] for j in range(p)) or any(lab[p + 1 + j] == lab[p + 2 + j] for j in range(q)):
                continue
            if not set(range(5)) <= set(lab):
                continue
            circ = Circle(twoarc_seq(lab[:p + 1], lab[p + 1:]), L, C, H, k=4)
            v = vec_from(circ); scale = max(1.0, np.abs(v).max())
            mg = F @ v / scale
            for o in range(len(ORB)):
                if o != 2:
                    worst[o] = min(worst[o], mg[ORB[o]].min())
    print('T1 two-arc circular-decomposable: min normalised margin per orbit', np.array2string(worst, precision=2), flush=True)
    # T2
    kept = 0; tried = 0; worst = np.full(len(ORB), np.inf); worst_slack_rejected = []
    for s in range(nsamp):
        p = int(rng.integers(2, 4)); q = int(rng.integers(2, 4))
        L, C, H = twoarc_circ(p, q, rng)
        n = p + q
        if rng.random() < 0.5:   # generic small noise
            eps = rng.choice([1e-3, 1e-2, 0.1])
            L2 = L + eps * np.triu(rng.normal(size=(n, n)), 1); L2 = np.triu(L2, 1); L2 = L2 + L2.T
            C2 = [np.maximum(0, C[0] + eps * rng.normal(size=n)), np.maximum(0, C[1] + eps * rng.normal(size=n))]
            h = max(0.0, H[0, 1] + eps * rng.normal()); H2 = np.array([[0, h], [h, 0]])
            L2 = np.maximum(L2, 0)
        else:                    # one-sided inflation of raw chord lengths (breaks raw Kalmanson)
            msk = rng.random((n, n)) < 0.4
            L2 = L + np.triu(msk * rng.exponential(size=(n, n)), 1); L2 = np.triu(L2, 1); L2 = L2 + L2.T
            C2 = [C[a] + (rng.random(n) < 0.4) * rng.exponential(size=n) for a in range(2)]
            h = H[0, 1] + (rng.random() < 0.4) * rng.exponential(); H2 = np.array([[0, h], [h, 0]])
        tried += 1
        sl, pur = physical_slack(lambda labs: twoarc_seq(labs[:p + 1], labs[p + 1:]), p + q + 2, L2, C2, H2, k=3)
        if sl < -1e-12:
            continue
        kept += 1
        for _ in range(8):
            lab = [int(x) for x in rng.integers(0, 6, size=p + q + 2)]
            if any(lab[j] == lab[j + 1] for j in range(p)) or any(lab[p + 1 + j] == lab[p + 2 + j] for j in range(q)):
                continue
            if not set(range(5)) <= set(lab):
                continue
            circ = Circle(twoarc_seq(lab[:p + 1], lab[p + 1:]), L2, C2, H2, k=3)
            v = vec_from(circ); scale = max(1.0, np.abs(v).max())
            mg = F @ v / scale
            for o in range(len(ORB)):
                if o != 2:
                    worst[o] = min(worst[o], mg[ORB[o]].min())
    print(f'T2 perturbed two-arc data passing all physical SSA/WM: kept {kept}/{tried}; min margin per orbit',
          np.array2string(worst, precision=2), flush=True)
