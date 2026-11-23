"""Laurent test at d13 = 0 of the residue-sum function Phi(d13) = sum_n Res_{d12=-n}[...] (WITHOUT Gamma(d13)^2):
   a0 = (1/2 pi i) oint Phi dd13 (non-zero iff Phi has a pole at 0, which would add a log^2 y term), and Phi(0) itself.
   Usage: python3 rh_t1_laurent.py <k> <dphi> <class> <amp> [N]"""
import sys
import numpy as np
import rh_common as R
import rh_t1_collinear as T

k, dphi, cl, nm = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3], sys.argv[4]
N = int(sys.argv[5]) if len(sys.argv) > 5 else (160 if cl in ('1243', '1324') else 12)
a, b, p2, p4 = T.cfg(dphi)
Sab = R.setup(k, a, b, p2, p4, classes=(cl,)); Sba = R.setup(k, b, a, p2, p4, classes=(cl,))
f = R.poly_D(int(nm[1:]))
out = []
for rad in (0.1, 0.2):
    th = 2*np.pi*np.arange(48)/48
    a0 = 0; mx = 0
    for t in th:
        d = rad*np.exp(1j*t)
        v = 0
        for o, S in (('ab', Sab), ('ba', Sba)):
            sg = 1.0 if o == 'ab' else (-1.0)**k
            v += sg*np.sum(R.residues(S, f, d, N)[cl])
        a0 += v*d/48; mx = max(mx, abs(v))
    out.append((rad, a0, mx))
    print(f'k {k} dphi {dphi} {cl} {nm}: radius {rad}: (1/2pi i) oint Phi = {a0:.3e}  (max|Phi| on circle {mx:.3e}; ratio {abs(a0)/mx:.1e})', flush=True)
