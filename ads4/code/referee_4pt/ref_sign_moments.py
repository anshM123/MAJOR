"""Referee M2 (round 2): K_odd/K_even from the exact cached modes with the referee's own weights:
int dth sin(th_out-th_in) G = (G_(-1,+1) - G_(+1,-1))/(2i)·(-1)?  derived:  G(th) = sum_q G_q e^{i q.th};
int sin(th3 - th1) e^{i(q1 th1 + q3 th3)} -> nonzero for (q1,q3) = (1,-1): 1/(2i)*... computed numerically here on a fine grid to avoid sign slips."""
import sys, os, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'fourpt')); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
from tmoments4 import TensorMoments4
from t5_analysis_l1 import hel_vectors
from cases import info
T = np.linspace(0, 2*np.pi, 64, endpoint=False)
def weights(f):
    # w(q1,q3) = mean over th1,th3 of f(th3 - th1) e^{i(q1 th1 + q3 th3)}
    A, B = np.meshgrid(T, T, indexing='ij')
    return {(q1, q3): np.mean(f(B - A)*np.exp(1j*(q1*A + q3*B))) for q1 in (-1, 0, 1) for q3 in (-1, 0, 1)}
WE, WO = weights(np.cos), weights(np.sin)
def K(mom, w, c, k, nu):
    tot = 0j
    for qL, wl in w.items():
        for qR, wr in w.items():
            if abs(wl*wr) < 1e-12: continue
            q = (qL[0], qL[1], qR[0], qR[1])
            for a in range(2):
                for b in range(2):
                    if c[a]*c[b] == 0: continue
                    v, sh = mom.moment(k, nu, (a, b), q)
                    tot += wl*wr*c[a]*c[b]*complex(v)*np.exp(sh - 700)
    return tot*(-1)**mom.J
for f in sorted(glob.glob(os.path.join(HERE, '..', 'fourpt', 'cache', 'modes4', '*.npz'))):
    mom = TensorMoments4.load(f); m, beta, _ = info(mom.Delta, mom.J)
    if float(m) < 60: continue
    o, s = hel_vectors(mom.Delta, mom.J)
    for nu in (0.0, 3.0):
        line = f"beta={float(beta):.2f} m={float(m):.1f} nu={nu}:"
        for name, c in (('same', s), ('opp', o), ('v1', (1.0, 0.0)), ('v2', (0.0, 1.0))):
            ke, ko = K(mom, WE, c, 6, nu), K(mom, WO, c, 6, nu)
            line += f"  {name}: Kodd/Keven = {ko.real/ke.real:+.5f} (Im/Re {abs(ko.imag/ko.real):.0e}, Keven sign {np.sign(ke.real):+.0f})"
        print(line, flush=True)
