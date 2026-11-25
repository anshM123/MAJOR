"""LF4: exact-de-shift cubic coefficient K(L) (reviewer formula, from Q6; checked to reproduce K = 0.99996752314 for NF02 at L = 998.5)
K(L) = -2mu26 + 245mu24/(2L^2) - 9619mu22/(8L^4) + 51575mu20/(32L^6) + 2mu34 - 35mu32/L^2 + 569mu30/(8L^4),  mu_kj = int_0^1 p^j psi_k,
L = Lam' = a Lam (support edge in nu units).  tau_K = tau_{K=1}/sqrt(K) (k(0) normalised to 1)."""
import json, glob, numpy as np
from fractions import Fraction as Fr
from ln_common import nf02_exact, mom_exact, guard
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'robust'))
from nf_lp import basis
def K(m2, m3, L):
    return (-2*m2(6) + 245*m2(4)/(2*L**2) - 9619*m2(2)/(8*L**4) + 51575*m2(0)/(32*L**6)
            + 2*m3(4) - 35*m3(2)/L**2 + 569*m3(0)/(8*L**4))
p2, p3 = nf02_exact(); k0 = -2*mom_exact(p2, 6) + 2*mom_exact(p3, 4)
out = []
for L in (Fr(9985, 10), Fr(9985, 100), Fr(49925, 1000), Fr(749, 10), Fr(2496, 100)):
    v = K(lambda j: mom_exact(p2, j)/k0, lambda j: mom_exact(p3, j)/k0, L)
    out.append('NF02 dilated, L = %s: K = %.11f' % (float(L), float(v)))
for f in sorted(glob.glob('cache/lp_a*_dl0.20_m1.json')):
    I = json.load(open(f)); a = I['a']; v = np.array(I['v']); n = I['d'] + 1; b = basis(I['d'], I['r'])
    Lb = max(len(x) for x in b); c2 = np.zeros(Lb); c3 = np.zeros(Lb)
    for i in range(n): c2[:len(b[i])] += v[i]*b[i]; c3[:len(b[i])] += v[n+i]*b[i]
    mm = lambda c: (lambda j: sum(ci/(i+j+1) for i, ci in enumerate(c)))
    if a == 1.0: continue
    Kv = K(mm(c2), mm(c3), a*998.5)
    out.append('LP a=%.3f delta=0.2: K(L=%.2f) = %.6f, tau(K=1) = %.4e, tau(K) = %.4e' % (a, a*998.5, Kv, I['tau'], I['tau']/np.sqrt(Kv)))
txt = '\n'.join(out); print(txt); open(guard('cache/ln7_K.txt'), 'w').write(txt + '\n')
