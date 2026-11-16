"""Compute and cache (cache/modes4/) the Fourier modes of the exact <JJJJ> partial waves G^{ab} on the moment grid.
Usage: python3 t4_compute_l1.py beta [m1 m2 ...]   (cases as in e2e/cases.py: J integer, (Delta, J) at fixed beta)."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'e2e'))
from cases import cases_fixed_beta          # noqa: E402
from tmoments4 import get_moments
beta = float(sys.argv[1])
ms = [float(x) for x in sys.argv[2:]] or [10, 20, 40, 70, 100]
for (D, J) in cases_fixed_beta(beta, ms):
    t = time.time()
    get_moments(1, D, J)
    print(f"beta={beta} Delta={D:.5f} J={J}: {time.time()-t:.0f} s", flush=True)
