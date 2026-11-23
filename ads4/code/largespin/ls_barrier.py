"""
largespin/ls_barrier.py (LARGESPIN_ANALYTIC.md §3): the impact-parameter (H2-density) form of the per-state action and the barrier
tolerance.

For an image of the k-th piece written as  I_k(nu) = int P_{-1/2+i nu}(cosh t) sigma_k(dt)  (on nu in [0, Lam]), Fubini gives exactly
    int_0^1 psi_k(p) I_k(Lam p) dp = int OmLR_k(t) sigma_k(dt).
A beyond-hybrid correction whose density sigma_d is carried by t >= t0 therefore moves the k-th piece by at most
    sup_{t >= t0} |OmLR_k(t)| * ||sigma_d||_TV .
Normalisation: the leading-Regge k = 2 image of a state is 2 P(cosh beta) (weight 2 at beta), the k = 3 image u sigma_h P (weight u).
Per-state barrier tolerances are computed in ls_struct.py:
    tau_k(state; eps) = floor / ( w_k * sup_{t >= (1-eps) beta} |OmLR_k(t)| ),   floor = N_flat/5 + min(0, curvature)  [certified + P/NE]
with w_2 = 2, w_3 = u.  Output: cache/barrier_<Dg>.txt, cache/omlr_<Dg>.npz.
Usage: python3 ls_barrier.py [Dg]
"""
import os, sys, json, time
sys.dont_write_bytecode = True
import numpy as np
import ls_core as C

Dg = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
L = C.lam(Dg)
fn = os.path.join(C.CACHE, 'omlr_%d.npz' % Dg)
if os.path.exists(fn):
    d = np.load(fn); ts, O2, O3, Of2 = d['ts'], d['O2'], d['O3'], d['Of2']
else:
    ts = np.concatenate([np.geomspace(2e-3, 0.02, 40, endpoint=False), np.geomspace(0.02, 4.0, 260)])
    t0 = time.time()
    O2 = C.OmLR(2, ts, Dg); O3 = C.OmLR(3, ts, Dg); Of2 = C.OmFlat(2, L*2*np.sinh(ts/2), Dg)
    np.savez(fn, ts=ts, O2=O2, O3=O3, Of2=Of2)
    print('OmLR table %.0f s' % (time.time() - t0))


def supabove(O, t0):
    m = ts >= t0
    return np.abs(O[m]).max() if m.any() else 0.0


out = []
w = out.append
w('Delta_g = %d, Lam = %.4f\n' % (Dg, L))
w('OmLR_2(t) (k = 2 leading-Regge action at impact parameter t, per unit ring weight); flat Regge analogue at b = 2 Lam sinh(t/2):')
w('%9s %9s %13s %13s %13s %10s' % ('t', 'bhat', 'OmLR_2', 'OmFlat_2', 'OmLR_3', 'b^3 OmLR2'))
for t in [0.003, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0]:
    i = int(np.argmin(abs(ts - t)))
    b = 2*L*np.sinh(ts[i]/2)
    w('%9.4f %9.2f %13.5e %13.5e %13.5e %10.1f' % (ts[i], b, O2[i], Of2[i], O3[i], b**3*O2[i]))
neg = ts[O2 <= 0]
w('\nOmLR_2 sign: min over t-grid = %.4e; negative at %d of %d grid points %s' % (O2.min(), len(neg), len(ts), (neg[:5] if len(neg) else '')))
d2 = np.diff(np.abs(O2))
w('|OmLR_2| non-increasing on the grid for t >= %.4f (first increase at t = %s)' %
  (ts[np.where(d2 > 0)[0][-1] + 1] if (d2 > 0).any() else ts[0], ts[np.where(d2 > 0)[0][:3] + 1] if (d2 > 0).any() else '-'))

# per-state tolerances: see ls_struct.py
txt = '\n'.join(out)
print(txt)
open(os.path.join(C.CACHE, 'barrier_%d.txt' % Dg), 'w').write(txt + '\n')
