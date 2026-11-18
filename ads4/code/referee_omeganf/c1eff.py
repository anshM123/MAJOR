"""Effective C1 implied by the skeleton scan's own 1/m^2 remainder for Omega^{NF,SC}: |S|/(u Ns), (E)/(u Ns), by bhat bin, Dg; also restricted
to fully controlled states (U == 0) -- the 'S' is the signed per-node skeleton remainder on controlled nodes."""
import pickle, os, numpy as np
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'omegau', 'cache')
for Dg in (20, 50, 100, 300, 1000):
    for Df in ('Df3', 'Df2'):
        R = pickle.load(open(os.path.join(C, f'recs_nf4_{Dg}_{Df}.pkl'), 'rb'))
        row = []
        for lo, hi in ((0, 10), (10, 20), (20, 41)):
            for ctl in (False, True):
                v = [r for h in R for r in R[h] if lo <= (2*r['J'] + 1)/r['mfl'] < hi and (r['U'] == 0 or not ctl)]
                if not v:
                    row.append(f"[{lo},{hi}){'c' if ctl else ''}: -"); continue
                s = [abs(r['S'])*r['mfl']**2/r['Ns'] for r in v]; e = [r['E']*r['mfl']**2/r['Ns'] for r in v]
                i = int(np.argmax(s))
                row.append(f"[{lo},{hi}){'c' if ctl else ' '} n={len(v):4d} max|S|/uNs={max(s):8.3g} (m={v[i]['mfl']:.2f},J={v[i]['J']}) med={np.median(s):.2g} maxE/uNs={max(e):8.3g}")
        print(f"Dg={Dg} {Df}:\n   " + "\n   ".join(row))
