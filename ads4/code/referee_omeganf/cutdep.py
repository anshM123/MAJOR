"""Skeleton threshold n(kappa_req, delta) (nf7 definitions) restricted to bhat <= B and to 10 < bhat <= B: cut dependence of the combined n."""
import pickle, os, numpy as np
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'omegau', 'cache')
NG = np.round(np.concatenate([np.arange(1.0, 1.3, 0.01), np.arange(1.3, 3.0, 0.05), [3.0]]), 4)
def marg(r, kq, dl): return r['Fh'] - r['E'] - kq*r['U'] - dl*r['Ns']
def n_of(L, kq, dl):
    for n0 in NG:
        if all(marg(r, kq, dl) > 0 for r in L if r['mfl'] >= n0 - 1e-9): return n0
    return np.inf
for Dg in (50, 100, 1000):
    R = [r for Df in ('Df3', 'Df2') for h in ('++', '+-') for r in pickle.load(open(os.path.join(C, f'recs_nf4_{Dg}_{Df}.pkl'), 'rb'))[h]]
    bh = lambda r: (2*r['J'] + 1)/r['mfl']
    row = []
    for lo, B in ((0, 10), (0, 11), (0, 12), (0, 15), (10, 12), (10, 15), (10, 20)):
        L = [r for r in R if lo < bh(r) <= B or (lo == 0 and bh(r) <= B)]
        row.append(f"({lo},{B}]: {n_of(L, 0.3, 0.03):.2f}/{n_of(L, 1.0, 0.0):.2f}")
    print(f"Dg={Dg}: n(0.3,0.03)/n(1,0): " + "  ".join(row))
