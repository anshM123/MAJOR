"""Adversarial hill-climb: can 'free-arc matching' data obey every physical
SSA / weak-monotonicity instance (all unions of asymptotic arcs) and still
violate a balanced 5-party facet?  Moves are random perturbations of the raw
chord data; a move is accepted only if all physical SSA/WM slacks stay >= 0
and the target (min normalised facet margin over a random set of labellings)
does not increase.
Run: python3 -I abstract_hillclimb.py geometry(one|two) nchains nsteps seed
"""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from abstract_axis import F, ORB, SUBSETS, circ_decomp
from abstract_twoarc import twoarc_seq, onearc_seq, twoarc_circ
from circle_engine import Circle

BAL = np.array([o != 2 for o in range(len(ORB))])
FB = F[np.concatenate([ORB[o] for o in range(len(ORB)) if o != 2])]


def all_unions_S(builder, nar, L, C, H, k):
    circ = Circle(builder(list(range(nar))), L, C, H, k=k)
    S = {}
    for r in range(0, nar + 1):
        for U in itertools.combinations(range(nar), r):
            U = frozenset(U); S[U] = circ.S(U) if U else 0.0
    return S


def phys_slack(S):
    keys = list(S); w = np.inf
    for U in keys:
        for V in keys:
            w = min(w, S[U] + S[V] - S[U | V] - S[U & V], S[U] + S[V] - S[U - V] - S[V - U])
    return w


def party_vec(S, lab):
    out = np.empty(31)
    for i, P in enumerate(SUBSETS):
        U = frozenset(j for j, l in enumerate(lab) if l in P)
        out[i] = S[U]
    return out


def target(S, labs):
    best = np.inf
    for lab in labs:
        v = party_vec(S, lab); sc = max(1.0, np.abs(v).max())
        best = min(best, (FB @ v).min() / sc)
    return best


def rand_labs(nar, cuts, rng, nl=12):
    """cuts: set of j such that arcs j, j+1 are separated by an edge (must differ)."""
    labs = []
    while len(labs) < nl:
        lab = [int(x) for x in rng.integers(0, 6, size=nar)]
        if any(lab[j] == lab[j + 1] for j in cuts):
            continue
        if set(range(5)) <= set(lab):
            labs.append(lab)
    return labs


def main():
    geom = sys.argv[1]; nch = int(sys.argv[2]); nst = int(sys.argv[3]); seed = int(sys.argv[4])
    rng = np.random.default_rng(seed)
    overall = np.inf
    for ch in range(nch):
        if geom == 'one':
            m = int(rng.integers(5, 7)); nar = m + 1
            L, c = circ_decomp(m, rng, sparse=0.25); C = [c]; H = np.zeros((1, 1))
            builder = lambda labs, m=m: onearc_seq(labs, m)
            cuts = set(range(m)); k = 4
        else:
            p = int(rng.integers(2, 4)); q = int(rng.integers(2, 4)); nar = p + q + 2; m = p + q
            L, C, H = twoarc_circ(p, q, rng, sparse=0.25)
            builder = lambda labs, p=p: twoarc_seq(labs[:p + 1], labs[p + 1:])
            cuts = set(range(p)) | set(p + 1 + j for j in range(q)); k = 3
        labs = rand_labs(nar, cuts, rng)
        S = all_unions_S(builder, nar, L, C, H, k)
        cur = target(S, labs); acc = 0
        for st in range(nst):
            sc = rng.choice([0.3, 0.1, 0.03, 0.01])
            L2 = L.copy(); C2 = [x.copy() for x in C]; H2 = H.copy()
            kind = rng.integers(0, 3)
            if kind == 0:
                i, j = rng.choice(m, 2, replace=False); d = sc * rng.normal()
                L2[i, j] = L2[j, i] = max(0, L2[i, j] + d)
            elif kind == 1:
                a = rng.integers(0, len(C)); i = rng.integers(0, m)
                C2[a][i] = max(0, C2[a][i] + sc * rng.normal())
            else:
                if H.shape[0] > 1:
                    h = max(0, H2[0, 1] + sc * rng.normal()); H2[0, 1] = H2[1, 0] = h
                else:
                    continue
            S2 = all_unions_S(builder, nar, L2, C2, H2, k)
            if phys_slack(S2) < -1e-12:
                continue
            t2 = target(S2, labs)
            if t2 <= cur + 1e-15:
                L, C, H, S, cur = L2, C2, H2, S2, t2; acc += 1
        overall = min(overall, cur)
        print(f'chain {ch} geom={geom} m={m} accepted={acc} final min normalised balanced margin={cur:+.3e}', flush=True)
    print('OVERALL min margin', overall)


if __name__ == '__main__':
    main()
