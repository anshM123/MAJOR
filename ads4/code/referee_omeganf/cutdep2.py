import pickle, os, numpy as np
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cutdep.py')).read().split('for Dg in')[0])
for Dg in (100, 1000):
    R = [r for Df in ('Df3', 'Df2') for h in ('++', '+-') for r in pickle.load(open(os.path.join(C, f'recs_nf4_{Dg}_{Df}.pkl'), 'rb'))[h]]
    bh = lambda r: (2*r['J'] + 1)/r['mfl']
    print(Dg, " ".join(f"({a},{b}]: {n_of([r for r in R if a < bh(r) <= b], 0.3, 0.03):.2f}/{n_of([r for r in R if a < bh(r) <= b], 1.0, 0.0):.2f}" for a, b in ((15, 20), (20, 30), (30, 41))))
