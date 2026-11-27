"""Vectorised exact (mod p) version of n8_five_point_count.py for larger D (18, 20).
Selects a spanning subset of monomials by rank on random points, then counts S5-symmetric f = P/<12>^8 (local N=8 MHV
five-point superamplitudes) at mass dimension D."""
import sys, itertools, time
import numpy as np, flint
import n8_five_point_count as base
PR = base.PR
rng = np.random.default_rng(11)

def kin_batch(n):
    out = []
    for _ in range(n):
        out.append(base.kin_mod())
    return out

def eval_matrix(mons, kins, perm=None):
    A = np.array([m[0] for m in mons], dtype=np.int64); B = np.array([m[1] for m in mons], dtype=np.int64)
    maxe = int(max(A.max(), B.max()))
    rows = []
    for lam, lt in kins:
        if perm is not None:
            lam = [lam[i] for i in perm]; lt = [lt[i] for i in perm]
        ang, sq = base.br_mod(lam, lt)
        # power tables
        val = np.ones(len(mons), dtype=np.int64)
        for k in range(10):
            pa = np.ones(maxe + 1, dtype=np.int64); ps = np.ones(maxe + 1, dtype=np.int64)
            for e in range(1, maxe + 1):
                pa[e] = pa[e - 1]*ang[k] % PR; ps[e] = ps[e - 1]*sq[k] % PR
            val = val*pa[A[:, k]] % PR
            val = val*ps[B[:, k]] % PR
        inv8 = pow(pow(ang[0], 8, PR), PR - 2, PR)
        rows.append(val*inv8 % PR)
    return np.array(rows, dtype=np.int64)

def rank_mod(M):
    return flint.nmod_mat(M.shape[0], M.shape[1], [int(x) for x in M.ravel()], PR).rank()

def count(D, probe=400):
    t0 = time.time()
    mons = base.monomials(D)
    print('D=%d: %d monomials (%.0fs)' % (D, len(mons), time.time() - t0)); sys.stdout.flush()
    kins = kin_batch(probe)
    E = eval_matrix(mons, kins)
    # greedy column selection via rank of random linear combos: use flint rref to pick pivot columns
    M = flint.nmod_mat(E.shape[0], E.shape[1], [int(x) for x in E.ravel()], PR)
    R, rk = M.rref()
    piv = []
    c = 0
    for r in range(rk):
        while int(R[r, c]) == 0: c += 1
        piv.append(c); c += 1
    print('  function-space dimension %d (probe %d points) (%.0fs)' % (rk, probe, time.time() - t0)); sys.stdout.flush()
    assert rk < probe - 20, 'increase probe'
    sub = [mons[i] for i in piv]
    kins2 = kin_batch(rk + 40)
    E2 = eval_matrix(sub, kins2)
    C = np.vstack([(eval_matrix(sub, kins2, perm=p) - E2) % PR for p in [(1, 0, 2, 3, 4), (1, 2, 3, 4, 0)]])
    rE, rC = rank_mod(E2), rank_mod(C)
    return len(mons), rE, rE - rC

if __name__ == '__main__':
    for D in [int(x) for x in sys.argv[1:]]:
        nm, r, ns = count(D)
        print('[exact mod p] D=%2d: monomials %6d, independent (--+++) polynomials %5d, S5-symmetric local N=8 superamplitudes (MHV sector): %d' % (D, nm, r, ns))
        sys.stdout.flush()
