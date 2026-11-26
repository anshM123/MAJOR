"""
kernel3.py: per-cell bound of |K_J(nu)| on the real axis (ANALYTIC_SLIVER3 §4), = rt2_kernel.K_cell (ANALYTIC_SLIVER2 §4: Theorem E',
Lemma X, Lemma F base powers, helicity for b = 4, Lemma B chain, REC tail with general growth G) with two changes:
  * the twist range is restricted to 998 <= a <= A_STAR = 4500: every state with a >= a*(L) is covered by the high-twist cell
    certificate (ANGULAR_ROUTE §73: a* = 4500 for beta_c in [3.04, 8.24], i.e. L in about [2.3, 7.3]; an imported input);
  * a finer, adaptive partition of the a-range (unit boxes near the turning point, then relative width `ratio`).
All bounds are those of rt2_kernel (imported read-only); only the cover of the a-range differs, and a cover is all the max needs.
"""
import math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'analytic2')); sys.path.insert(0, os.path.join(HERE, '..', 'analytic'))
import common2
import rt2_kernel as RK
A_STAR = 4500.0


def a_ranges(T1, ratio=1.03, a_max=A_STAR):
    out = []; a = 998.0
    while a < a_max:
        w = 1.0 if a < T1 + 30 else max(1.0, a*(ratio - 1))
        b = min(a + w, a_max)
        out.append((a, b)); a = b
    return out


def K_cell(T0, T1, La, Lb, ratio=1.03, use_frames=True, jmax=110, J0x=1500, verbose=False):
    return RK.K_cell(T0, T1, La, Lb, verbose=verbose, use_frames=use_frames, ranges=a_ranges(T1, ratio), jmax=jmax, J0x=J0x)


if __name__ == '__main__':
    import time
    T0, La, Lb = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    for ratio in (1.15, 1.03):
        t0 = time.time()
        print(ratio, K_cell(T0, T0 + 1.1, La, Lb, ratio=ratio), f'{time.time()-t0:.0f}s', flush=True)
