"""which top-Regge degree J is borderline (first non-zero symmetrised term-wise value) in each Wick class, k = 2 and 4 (P1 angles)."""
import sys, time, pickle
import numpy as np
import rh_common as R

k = int(sys.argv[1])
a, b, p2, p4 = R.P.config('P1')
t0 = time.time()
Sab = R.setup(k, a, b, p2, p4); Sba = R.setup(k, b, a, p2, p4)
amps = {f'D{J}': R.poly_D(J) for J in range(2, 9)}
out, acc = R.values(Sab, Sba, k, amps, N=160)
lines = [f'k = {k}, P1, symmetrised term-wise values of D^J/4 (D = d12 - d14), per class; N = 160 + tail [{time.time()-t0:.0f}s]']
for J in range(2, 9):
    for c in R.C.CLS:
        v, tl, pp = out[(f'D{J}', c)]
        lines.append(f'  J = {J}  {c:5s}  {v.real:+.6e} {v.imag:+.1e}i   tail {np.real(tl):+.1e} (p {pp:.2f})')
open(f'{R.CACHE}/degree_k{k}.txt', 'w').write('\n'.join(lines) + '\n')
pickle.dump(acc, open(f'{R.CACHE}/degree_k{k}.pkl', 'wb'))
print('\n'.join(lines))
